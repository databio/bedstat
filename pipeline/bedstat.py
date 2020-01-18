#!/usr/bin/env python
from argparse import ArgumentParser
import pypiper, os, sys
from elasticsearch import Elasticsearch
from elasticsearch import helpers
from elasticsearch.serializer import JSONSerializer
import json

parser = ArgumentParser(
    description="A pipeline to read a file in BED format and produce metadata in JSON format.")

parser.add_argument('--bedfile', help='full path to bed file to process', required=True)
# parser.add_argument('--outfolder', default=output_dir, help='folder to put images and json files in')

parser = pypiper.add_pypiper_args(parser, args=["genome"], groups=["pypiper", "common", "looper", "ngs"],
                                  required=['bedfile', 'genome'])

args = parser.parse_args()

outfolder = args.output_parent

# get the --fileid argument for the R script
split_path = os.path.split(args.bedfile)
bedfile_portion = split_path[1]

fileid = os.path.splitext(bedfile_portion)[0]

# get the output folder argument for the R script
outfolder = os.path.abspath(os.path.join(outfolder, fileid))

# try to create the directory and ignore failure if it already exists
# os.makedirs(outfolder, exist_ok=True)

pm = pypiper.PipelineManager(name="bedstat-pipeline", outfolder=outfolder, args=args)

rscript_path = os.path.join(os.path.pardir(os.path.dirname(os.path.abspath(__file__))), "tools", "regionstat.R")

cmd_vars = dict(rscript=rscript_path, bed=args.bedfile, id=fileid, out=outfolder, genome=args.genome_assembly)
command = "Rscript {rscript} --bedfile={bed} --fileid={id} --outputfolder={out} --genome={genome}".format(cmd_vars)

target = os.path.abspath(os.path.join(outfolder, bedfile_portion))

pm.run(command, target)

# now get the resulting json file and load it into elasticsearch
# it the file exists, of course
if os.path.splitext(bedfile_portion)[1] != '':
    # open connection to elastic
    try:
        es = Elasticsearch([{'host': 'localhost'}])
        json_file = os.path.splitext(bedfile_portion)[0] + ".json"
        json_file_path = os.path.abspath(os.path.join(outfolder, json_file))
        with open(json_file_path, 'r', encoding='utf-8') as f:
            data = json.loads(f.read())
        es.index(index="bedstat_bedfiles", doc_type='doc', body=data)
    except Exception as e:
        print(e)

pm.stop_pipeline()

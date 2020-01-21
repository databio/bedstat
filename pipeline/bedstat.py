#!/usr/bin/env python
from argparse import ArgumentParser
from elasticsearch import Elasticsearch
import json
import gzip
import shutil
import yaml
import os
import sys

import pypiper
from bbconf.const import *

parser = ArgumentParser(
    description="A pipeline to read a file in BED format and produce metadata in JSON format.")

parser.add_argument('--bedfile', help='full path to bed file to process', required=True)
parser.add_argument('--nodbcommit', help='do not try to commit json outout to back end database', action='store_true')
parser.add_argument('--dbhost',
                    help='use database host address/name to connect to',
                    required='--nodbcommit' not in sys.argv,
                    default='localhost')
# parser argument for yaml - to pass on more sample metadata into the database
parser.add_argument("-y", "--sample-yaml",
                    dest="sample_config",
                    help="Yaml config file with sample attributes.",
                    type=str)
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

# get the sample line from the yaml config file
y = yaml.safe_load(open(args.sample_config, "r"))

pm = pypiper.PipelineManager(name="bedstat-pipeline", outfolder=outfolder, args=args)

rscript_path = os.path.join(os.path.pardir(os.path.dirname(os.path.abspath(__file__))), "tools", "regionstat.R")

cmd_vars = dict(rscript=rscript_path, bed=args.bedfile, id=fileid, out=outfolder, genome=args.genome_assembly)
command = "Rscript {rscript} --bedfile={bed} --fileid={id} --outputfolder={out} --genome={genome}".format(cmd_vars)

# gzip the original bed file and keep it with the pipeline results
try:
    dst_path = os.path.abspath(os.path.join(outfolder, bedfile_portion))
    with open(args.bedfile, 'rb') as f_in:
        with gzip.open(dst_path + '.gz', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    # create a symlink to make the .gz easier to find
    link_dest = os.path.abspath(os.path.join(outfolder, "raw_bedfile"))
    if not os.path.exists(link_dest):
        os.symlink(bedfile_portion + ".gz", link_dest)
except Exception as e:
    raise e

target = os.path.abspath(os.path.join(outfolder, bedfile_portion))

pm.run(command, target)

# now get the resulting json file and load it into elasticsearch
# it the file exists, of course
if not args.nodbcommit and os.path.splitext(bedfile_portion)[1] != '':
    # open connection to elastic
    try:
        es = Elasticsearch([{'host': args.dbhost}])
        json_file = os.path.splitext(bedfile_portion)[0] + ".json"
        json_file_path = os.path.abspath(os.path.join(outfolder, json_file))
        with open(json_file_path, 'r', encoding='utf-8') as f:
            data = json.loads(f.read())
            for key in SEARCH_TERMS:
                try:
                    data[key] = y[key]
                except KeyError:
                    pm.warning("Can't find key: {}".format(key))
                    pm.warning(y)
                # enrich the data from R with the data from the sample line itself
            pm.info("Inserting into database...")
            es.index(index=BED_INDEX, body=data)
        es.index(index=BED_INDEX, doc_type='doc', body=data)
    except Exception as e:
        print(e)

pm.stop_pipeline()

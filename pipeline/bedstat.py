#!/usr/bin/env python
from argparse import ArgumentParser
import pypiper, os, sys
from elasticsearch import Elasticsearch
from elasticsearch import helpers
from elasticsearch.serializer import JSONSerializer
import json, gzip, shutil, yaml

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

parser = pypiper.add_pypiper_args(parser, args=["genome"], groups=["pypiper", "common", "looper", "ngs"],
                                  required=['bedfile', 'genome'])

args = parser.parse_args()

outfolder = args.output_parent

# get the --bedfile argument for the R script
bfile = args.bedfile

# get the --fileid argument for the R script
split_path = os.path.split(bfile)

bedfile_portion = split_path[1]

fileid = os.path.splitext(bedfile_portion)[0]

# get the output folder argument for the R script
outfolder = os.path.abspath(os.path.join(outfolder, fileid))

# get the sample line from the yaml config file
y = yaml.load(open(args.sample_config, "r"))

pm = pypiper.PipelineManager(name="bedstat-pipeline", outfolder=outfolder, args=args)

command = "Rscript tools/regionstat.R --bedfile=%s --fileid=%s --outputfolder=%s --genome=%s" \
          % (bfile, fileid, outfolder, args.genome_assembly)

target = os.path.abspath(os.path.join(outfolder, bedfile_portion))

# gzip the original bed file and keep it with the pipeline results
try:
    dst_path = os.path.abspath(os.path.join(outfolder, bedfile_portion))
    with open(bfile, 'rb') as f_in:
        with gzip.open(dst_path + '.gz', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    # create a symlink to make the .gz easier to find
    if not os.path.exists("raw_bedfile"):
        os.symlink(bedfile_portion + ".gz", os.path.abspath(os.path.join(outfolder, "raw_bedfile")))
except Exception as e:
    raise e

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
            for key in ['cellType', 'cellTypeSubtype', 'antibody', 'mappingGenome', 'description', 'tissue', 'species']:
                data[key] = y[key]
            # enrich the data from R with the data from the sample line itself
        es.index(index="bedstat_bedfiles", body=data)
    except Exception as e:
        raise e

pm.stop_pipeline()

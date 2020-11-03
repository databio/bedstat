#!/usr/bin/env python3
"""
bedfile statistics generating pipeline
"""

__author__ = ["Michal Stolarczyk", "Ognen Duzlevski", "Jose Verdezoto"]
__email__ = "michal@virginia.edu"
__version__ = "0.0.2-dev"

from argparse import ArgumentParser
from hashlib import md5
from psycopg2.extras import Json
import json
import yaml
import os
import warnings

from bbconf.const import *
import pypiper
import bbconf

parser = ArgumentParser(
    description="A pipeline to read a file in BED format and produce metadata "
                "in JSON format.")

parser.add_argument(
    '--bedfile', help='a full path to bed file to process', required=True)
parser.add_argument(
    '--open-signal-matrix', type=str, required=False, default=None,
    help='a full path to the openSignalMatrix required for the tissue '
         'specificity plots')
parser.add_argument(
    "--bedbase-config", dest="bedbase_config", type=str, default=None,
    help="a path to the bedbase configuratiion file")
parser.add_argument(
    "-y", "--sample-yaml", dest="sample_yaml", type=str, required=False,
    help="a yaml config file with sample attributes to pass on more metadata "
         "into the database")
exclusive_group = parser.add_mutually_exclusive_group()
exclusive_group.add_argument(
    '--no-db-commit', action='store_true',
    help='whether the JSON commit to the database should be skipped')
exclusive_group.add_argument(
    '--just-db-commit', action='store_true',
    help='whether just to commit the JSON to the database')
parser = pypiper.add_pypiper_args(parser,
                                  groups=["pypiper", "common", "looper", "ngs"])

args = parser.parse_args()

bbc = bbconf.BedBaseConf(filepath=bbconf.get_bedbase_cfg(args.bedbase_config))

bed_digest = md5(open(args.bedfile, 'rb').read()).hexdigest()
bedfile_name = os.path.split(args.bedfile)[1]
# need to split twice since there are 2 exts
fileid = os.path.splitext(os.path.splitext(bedfile_name)[0])[0]
outfolder = os.path.abspath(os.path.join(
    bbc[CFG_PATH_KEY][CFG_BEDSTAT_OUTPUT_KEY], bed_digest))
json_file_path = os.path.abspath(os.path.join(outfolder, fileid + ".json"))
json_plots_file_path = os.path.abspath(os.path.join(outfolder,
                                                    fileid + "_plots.json"))
bed_relpath = os.path.relpath(
    args.bedfile, os.path.abspath(bbc[CFG_PATH_KEY][CFG_BEDSTAT_OUTPUT_KEY]))

if not args.just_db_commit:
    pm = pypiper.PipelineManager(name="bedstat-pipeline", outfolder=outfolder,
                                 args=args)
    rscript_path = os.path.join(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__))), "tools", "regionstat.R")
    assert os.path.exists(rscript_path), \
        FileNotFoundError(f"'{rscript_path}' script not found")
    command = \
        f"Rscript {rscript_path} --bedfilePath={args.bedfile} " \
        f"--fileId={fileid} --openSignalMatrix={args.open_signal_matrix} " \
        f"--outputFolder={outfolder} --genome={args.genome_assembly} " \
        f"--digest={bed_digest} --bedfileRelpath={bed_relpath}"
    pm.run(cmd=command, target=json_file_path)
    pm.stop_pipeline()

# now get the resulting json file and load it into Elasticsearch
# if the file exists, of course
if not args.no_db_commit:
    data = {}
    with open(json_file_path, 'r', encoding='utf-8') as f:
        data = json.loads(f.read())
    with open(json_plots_file_path, 'r', encoding='utf-8') as f_plots:
        plots = json.loads(f_plots.read())
    if args.sample_yaml:
        # get the sample-specific metadata from the sample yaml representation
        other = {}
        if os.path.exists(args.sample_yaml):
            y = yaml.safe_load(open(args.sample_yaml, "r"))
            for key in JSON_METADATA_VALUES:
                # keep just the metadata we care about
                try:
                    other[key] = y[key]
                except KeyError:
                    print(f"'{key}' metadata not available")
        else:
            warnings.warn(
                f"Specified sample_yaml path does not exist: {args.sample_yaml}")
    # enrich the data from R with the data from the sample line itself
    # unlist the data, since the output of regionstat.R is a dict of lists of
    # length 1 and force keys to lower to correspond with the
    # postgres column indentifiers
    data.update({JSON_OTHER_KEY: other})
    data = {k.lower(): v[0] if isinstance(v, list) else v for k, v in data.items()}
    data.update(dict(plots=Json(plots)))
    if not bbc.check_bedfiles_table_exists():
        bbc.create_bedfiles_table(columns=BED_COLUMNS)
    bbc.insert_bedfile_data(values=data)

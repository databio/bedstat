#!/usr/bin/env python3
"""
bedfile statistics generating pipeline
"""

__author__ = ["Michal Stolarczyk", "Ognen Duzlevski", "Jose Verdezoto"]
__email__ = "michal@virginia.edu"
__version__ = "0.0.3-dev"

from argparse import ArgumentParser
from hashlib import md5
import pandas as pd
import json
import yaml
import os
import sys
import warnings
import tempfile
import requests

import pypiper
import bbconf
import time

parser = ArgumentParser(
    description="A pipeline to read a file in BED format and produce metadata "
    "in JSON format."
)

parser.add_argument(
    "--bedfile", help="a full path to bed file to process", required=True
)
parser.add_argument(
    "--open-signal-matrix",
    type=str,
    required=False,
    default=None,
    help="a full path to the openSignalMatrix required for the tissue "
    "specificity plots",
)
parser.add_argument(
    "--bigbed",
    type=str,
    required=False,
    default=None,
    help="a full path to the bigbed files",
)
parser.add_argument(
    "--bedbase-config",
    dest="bedbase_config",
    type=str,
    default=None,
    help="a path to the bedbase configuratiion file",
)
parser.add_argument(
    "-y",
    "--sample-yaml",
    dest="sample_yaml",
    type=str,
    required=False,
    help="a yaml config file with sample attributes to pass on more metadata "
    "into the database",
)
exclusive_group = parser.add_mutually_exclusive_group()
exclusive_group.add_argument(
    "--no-db-commit",
    action="store_true",
    help="whether the JSON commit to the database should be skipped",
)
exclusive_group.add_argument(
    "--just-db-commit",
    action="store_true",
    help="whether just to commit the JSON to the database",
)
parser = pypiper.add_pypiper_args(parser, groups=["pypiper", "common", "looper", "ngs"])

args = parser.parse_args()

bbc = bbconf.BedBaseConf(config_path=args.bedbase_config, database_only=True)
bedstat_output_path = bbc.get_bedstat_output_path()

bed_digest = md5(open(args.bedfile, "rb").read()).hexdigest()
bedfile_name = os.path.split(args.bedfile)[1]
# need to split twice since there are 2 exts
fileid = os.path.splitext(os.path.splitext(bedfile_name)[0])[0]
outfolder = os.path.abspath(os.path.join(bedstat_output_path, bed_digest))
json_file_path = os.path.abspath(os.path.join(outfolder, fileid + ".json"))
json_plots_file_path = os.path.abspath(os.path.join(outfolder, fileid + "_plots.json"))
bed_relpath = os.path.relpath(
    args.bedfile, os.path.abspath(os.path.join(bedstat_output_path, os.pardir, os.pardir))
)
bigbed_relpath = os.path.relpath(
    os.path.join(args.bigbed, fileid + ".bigBed"),
    os.path.abspath(os.path.join(bedstat_output_path, os.pardir, os.pardir)),
)


def convert_unit(size_in_bytes):
    """ Convert the size from bytes to other units like KB, MB or GB"""
    if size_in_bytes < 1024:
        return str(size_in_bytes) + "bytes"
    elif size_in_bytes in range(1024, 1024*1024):
        return str(round(size_in_bytes/1024, 2)) + "KB"
    elif size_in_bytes in range(1024*1024, 1024*1024*1024):
        return str(round(size_in_bytes/(1024*1024))) + "MB"
    elif size_in_bytes >= 1024*1024*1024:
        return str(round(size_in_bytes/(1024*1024*1024))) + "GB"


def get_file_size(file_name):
   """ Get file in size in given unit like KB, MB or GB"""
   size = os.path.getsize(file_name)
   return convert_unit(size)


def main():
    if not args.just_db_commit:
        pm = pypiper.PipelineManager(
            name="bedstat-pipeline", outfolder=outfolder, args=args
        )

        # run Rscript
        rscript_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "tools",
            "regionstat.R",
        )
        assert os.path.exists(rscript_path), FileNotFoundError(
            f"'{rscript_path}' script not found"
        )
        command = (
            f"Rscript {rscript_path} --bedfilePath={args.bedfile} "
            f"--fileId={fileid} --openSignalMatrix={args.open_signal_matrix} "
            f"--outputFolder={outfolder} --genome={args.genome_assembly} "
            f"--digest={bed_digest}"
        )
        print(command)
        pm.run(cmd=command, target=json_file_path)
        pm.stop_pipeline()

    # now get the resulting json file and load it into Elasticsearch
    # if the file exists, of course
    if not args.no_db_commit:
        data = {}
        with open(json_file_path, "r", encoding="utf-8") as f:
            data = json.loads(f.read())
        with open(json_plots_file_path, "r", encoding="utf-8") as f_plots:
            plots = json.loads(f_plots.read())
        if args.sample_yaml:
            # get the sample-specific metadata from the sample yaml representation
            other = {}
            if os.path.exists(args.sample_yaml):
                y = yaml.safe_load(open(args.sample_yaml, "r"))
                data.update({"other": y})
        # unlist the data, since the output of regionstat.R is a dict of lists of
        # length 1 and force keys to lower to correspond with the
        # postgres column identifiers
        data = {k.lower(): v[0] if isinstance(v, list) else v for k, v in data.items()}
        data.update({"bedfile": {"path": bed_relpath, "size": get_file_size(args.bedfile), "title": "Path to the BED file"}})
        
        if os.path.exists(
            os.path.join(args.bigbed, fileid + ".bigBed")
            ) and not os.path.islink(os.path.join(args.bigbed, fileid + ".bigBed")):
            digest = requests.get(
                f"http://refgenomes.databio.org/genomes/genome_digest/{args.genome_assembly}"
            ).text.strip('""')

            data.update(
                {
                    "genome": {
                        "alias": args.genome_assembly,
                        "digest": digest,
                    }
                }
            )
            data.update(
                {
                    "bigbedfile": {
                        "path": bigbed_relpath, 
                        "size": get_file_size(os.path.join(args.bigbed, fileid + ".bigBed")),
                        "title": "Path to the big BED file",
                    }
                }
            )

        for plot in plots:
            plot_id = plot["name"]
            del plot["name"]
            data.update({plot_id: plot})
        bbc.bed.report(record_identifier=bed_digest, values=data)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Pipeline aborted.")
        sys.exit(1)

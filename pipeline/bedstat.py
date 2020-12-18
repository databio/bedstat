#!/usr/bin/env python3
"""
bedfile statistics generating pipeline
"""

__author__ = ["Michal Stolarczyk", "Ognen Duzlevski", "Jose Verdezoto"]
__email__ = "michal@virginia.edu"
__version__ = "0.0.2-dev"

from argparse import ArgumentParser
from hashlib import md5
import pandas as pd
import json
import yaml
import os
import warnings
import tempfile

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
    '--bigbed', type=str, required=False, default=None,
    help='a full path to the bigbed files')
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

bbc = bbconf.BedBaseConf(config_path=args.bedbase_config, database_only=True)
bedstat_output_path = bbc.get_bedstat_output_path()

bed_digest = md5(open(args.bedfile, 'rb').read()).hexdigest()
bedfile_name = os.path.split(args.bedfile)[1]
# need to split twice since there are 2 exts
fileid = os.path.splitext(os.path.splitext(bedfile_name)[0])[0]
outfolder = os.path.abspath(os.path.join(bedstat_output_path, bed_digest))
json_file_path = os.path.abspath(os.path.join(outfolder, fileid + ".json"))
json_plots_file_path = os.path.abspath(
    os.path.join(outfolder, fileid + "_plots.json"))
bed_relpath = os.path.relpath(
    args.bedfile, os.path.join(os.path.abspath(bedstat_output_path), bed_digest))
print ("check:", args.bigbed)
bigbed_relpath = os.path.relpath(
    os.path.join(args.bigbed, fileid + ".bigBed"), os.path.join(os.path.abspath(bedstat_output_path), bed_digest))

if not args.just_db_commit:
    pm = pypiper.PipelineManager(name="bedstat-pipeline", outfolder=outfolder,
                                 args=args)

    # # Produce bigBed (bigNarrowPeak) file from peak file 
    # bigNarrowPeak = os.path.join(outfolder, fileid + ".bigBed")
    # temp = tempfile.NamedTemporaryFile(dir=outfolder, delete=False)
    # print ("test bigbed saving path: ", bigNarrowPeak)
    # print ("test chrom.sizes path: ", args.chrom_size)
    # if not os.path.exists(bigNarrowPeak):
    #     df = pd.read_csv(args.bedfile, sep='\t', header=None,
    #                         names=("V1","V2","V3","V4","V5","V6",
    #                                 "V7","V8","V9","V10")).sort_values(by=["V1","V2"])
    #     df.to_csv(temp.name, sep='\t', header=False, index=False)
    #     pm.clean_add(temp.name)
    #     print ("BED: \n", df)
    #     as_file = os.path.join(outfolder, "bigNarrowPeak.as")
    #     cmd = ("echo 'table bigNarrowPeak\n" + 
    #             "\"BED6+4 Peaks of signal enrichment based on pooled, normalized (interpreted) data.\"\n" +
    #             "(\n" +
    #             "     string chrom;        \"Reference sequence chromosome or scaffold\"\n" +
    #             "     uint   chromStart;   \"Start position in chromosome\"\n" +
    #             "     uint   chromEnd;     \"End position in chromosome\"\n" +
    #             "     string name;         \"Name given to a region (preferably unique). Use . if no name is assigned\"\n" +
    #             "     uint   score;        \"Indicates how dark the peak will be displayed in the browser (0-1000) \"\n" +
    #             "     char[1]  strand;     \"+ or - or . for unknown\"\n" +
    #             "     float  signalValue;  \"Measurement of average enrichment for the region\"\n" +
    #             "     float  pValue;       \"Statistical significance of signal value (-log10). Set to -1 if not used.\"\n" +
    #             "     float  qValue;       \"Statistical significance with multiple-test correction applied (FDR -log10). Set to -1 if not used.\"\n" +
    #             "     int   peak;          \"Point-source called for this peak; 0-based offset from chromStart. Set to -1 if no point-source called.\"\n" +
    #             ")' > " + as_file)
    #     pm.run(cmd, as_file, clean=True)

    #     cmd = ("bedToBigBed -as=" + as_file + " -type=bed6+4 " +
    #             temp.name + " " + args.chrom_size + " " + bigNarrowPeak)
    #     pm.run(cmd, bigNarrowPeak, nofail=True)

    # run Rscript
    rscript_path = os.path.join(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__))), "tools", "regionstat.R")
    assert os.path.exists(rscript_path), \
        FileNotFoundError(f"'{rscript_path}' script not found")
    command = \
        f"Rscript {rscript_path} --bedfilePath={args.bedfile} " \
        f"--fileId={fileid} --openSignalMatrix={args.open_signal_matrix} " \
        f"--outputFolder={outfolder} --genome={args.genome_assembly} " \
        f"--digest={bed_digest}"
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
            data.update({"other": y})
    # unlist the data, since the output of regionstat.R is a dict of lists of
    # length 1 and force keys to lower to correspond with the
    # postgres column identifiers
    data = {k.lower(): v[0] if isinstance(v, list) else v for k, v in data.items()}
    data.update({"bedfile": {"path": bed_relpath, "title": "Path to the BED file"}})
    data.update({"bigbedfile": {"path": bigbed_relpath, "title": "Path to the big BED file"}})
    for plot in plots:
        plot_id = plot["name"]
        del plot["name"]
        data.update({plot_id: plot})
    bbc.bed.report(record_identifier=bed_digest, values=data)

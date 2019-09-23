#!/usr/bin/env python

from argparse import ArgumentParser
import pypiper, os, sys

parser = ArgumentParser(
    description="A pipeline to read a file in BED format and produce metadata in JSON format.")

parser.add_argument('--bedfile', help='full path to bed file to process')
parser.add_argument('--outfolder', default='output', help='folder to put images and json files in')

parser = pypiper.add_pypiper_args(parser, args=["genome"], groups=["pypiper", "common", "looper", "ngs"], required=['bedfile', 'genome'])

args = parser.parse_args()

# get all the arguments and prepare to call R script
if (args.bedfile is None):
    parser.print_help()
    sys.exit()

# get the --bedfile argument for the R script
bfile = args.bedfile

# get the --fileid argument for the R script
split_path = os.path.split(bfile)

bedfile_portion = split_path[1]
dot_separator_idx = bedfile_portion.find('.')
if (dot_separator_idx < 0):
    fileid = bedfile_portion
else:
    fileid = bedfile_portion[0:dot_separator_idx]

# get the output folder argument for the R script
outfolder = os.path.abspath(os.path.join(args.outfolder, fileid))

# try to create the directory and ignore failure if it already exists
#os.makedirs(outfolder, exist_ok=True)

pm = pypiper.PipelineManager(name="bedstat-pipeline", outfolder=outfolder, args=args)

command = "Rscript tools/regionstat.R --bedfile=%s --fileid=%s --outputfolder=%s --genome=%s" % (bfile, fileid, outfolder, args.genome_assembly)

target = os.path.abspath(os.path.join(outfolder, bedfile_portion))

pm.run(command, target)

pm.stop_pipeline()

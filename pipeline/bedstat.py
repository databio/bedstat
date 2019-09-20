#!/usr/bin/env python

from argparse import ArgumentParser
import pypiper, os

parser = ArgumentParser(
    description="A pipeline to read a file in BED format and produce metadata in JSON format.")

parser = pypiper.add_pypiper_args(parser, groups=["pypiper", "common", "looper"],
                                    args=["output-parent", "input"],
                                    required=['sample-name', 'output-parent', 'input'])

args = parser.parse_args()

#outfolder = "bed_to_json_output/"
outfolder = os.path.abspath(os.path.join(args.output_parent, args.sample_name))

pm = pypiper.PipelineManager(name="bed-to-json-pipeline", outfolder=outfolder, args=args)

raw_folder = os.path.join(outfolder, "raw/")

local_input_file = raw_folder + "/" + args.input

target = os.path.join(outfolder, "bedfile.out")
command = "Rscript regionstat.R"

pm.run(command, target)
pm.stop_pipeline()

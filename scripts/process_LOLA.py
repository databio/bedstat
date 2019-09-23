#!/usr/bin/env python

from argparse import ArgumentParser
import os, sys
from pathlib import Path

# for reading TSV files
import csv

# process index file from LOLA core DB
# expects a base directory for the collection
# and base_dir/regions to contain the bed files
# the index file is expected to be in tsv format
# however, the header structure is not fixed
def process_index_file(base_dir, idx_file):
    with open(idx_file) as tsvfile:
        reader = csv.DictReader(tsvfile, dialect='excel-tab')
        for row in reader:
            # do something with each row
            pass

def process_collection_file(base_dir, cfile):
    print("Processing collection file %s" % cfile)
    print("Base dir: %s" % base_dir)

if __name__ == '__main__':
    parser = ArgumentParser(description="Script to process LOLA database into PEP format.")

    parser.add_argument('--lola_loc', help='full path to base LOLA db (folder whose subfolders are hg38, hg19 etc.')
    parser.add_argument('--genome', help='genome to process, e.g. hg38', default='hg38')

    args = parser.parse_args()

    # get all the arguments and prepare to call R script
    if (args.lola_loc is None):
        parser.print_help()
        sys.exit()

    lola_base = args.lola_loc
    genome = args.genome

    # we are interested in index.txt and collection.txt files
    # each file has a header and is \t delimited, however, header structure is not fixed
    full_path = os.path.abspath(os.path.join(lola_base, genome))

    # search for subfolders and index.txt
    for filename in Path(full_path).glob('**/index.txt'):
        # break down the path from each index.txt file
        # to get its base directory relative to lola_loc
        # and see if collection.txt can be obtained as well
        try:
            path_parts = os.path.split(filename)
            if (len(path_parts) != 2): break
            collection_file = os.path.abspath(os.path.join(path_parts[0], 'collection.txt'))
            cpath = Path(collection_file)
            if (cpath.exists() and cpath.is_file()):
                # we also have a collection file
                process_index_file(path_parts[0], filename)
                process_collection_file(path_parts[0], cpath)
        except Exception as e:
            print("Error accessing file: {0}".format(e))


        

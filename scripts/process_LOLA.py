#!/usr/bin/env python

# purpose of this script is to run through LOLA Core DB and produce a CSV file as input for looper
# Author: Ognen Duzlevski

from argparse import ArgumentParser
import os, sys
from pathlib import Path

# for reading TSV files
import csv

class LOLAIndexFileError(Exception):
    pass

# process index file from LOLA core DB
# expects a base directory for the collection
# and base_dir/regions to contain the bed files
# the index file is expected to be in tsv format
# however, the header structure is not fixed
# if the file happens to be in csv format, we try to make it right
def process_index_file(base_dir, idx_file, genome):
    #print("%s" % idx_file)
    bed_path_base = os.path.abspath(os.path.join(base_dir, "regions"))
    with open(idx_file) as tsvfile, open(idx_file) as test:
        try:
            reader = csv.DictReader(test, dialect='excel-tab')
            n = reader.__next__()
            l = list(n)
            if len(n) == 1:
                # csv case
                if (l[0].find(',') != -1):
                    test.seek(0)
                    reader = csv.DictReader(test)
                else:
                    raise Exception("Not tsv and not csv format for %s. Aborting." % idx_file)
            else:
                # tsv case
                reader = csv.DictReader(tsvfile, dialect='excel-tab')
            for row in reader:
                s = ''
                if 'filename' in row:
                    # see if the file name makes sense - sometimes some filenames have , or . in them
                    # and are split inappropriately
                    s = os.path.abspath(os.path.join(bed_path_base, row['filename']))
                    f = Path(s)
                    if not (f.exists() and f.is_file()):
                        # skip file
                        raise LOLAIndexFileError
                else:
                    raise LOLAIndexFileError
                if 'cellType' in row:
                    s = s + "," + row['cellType']
                else:
                    s = s + ',' + "unknown"
                s = s + ',bedstat' + ',' + genome
                print(s)
        except LOLAIndexFileError as l:
            sys.stderr.write("Skipping line %s in index file %s\n" % (row, idx_file))
            pass
        except Exception as e:
            print("Exception reading %s file" % idx_file)
            raise

def process_collection_file(base_dir, cfile):
    #print("Processing collection file %s" % cfile)
    #print("Base dir: %s" % base_dir)
    pass

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

    # what we are interested in for now from index.txt files
    print("sample_name,celltype,protocol,genome")
    # search for subfolders and index.txt
    for filename in Path(full_path).glob('**/index.txt'):
        # break down the path from each index.txt file
        # to get its base directory relative to lola_loc
        # and see if collection.txt can be obtained as well
        try:
            path_parts = os.path.split(filename)
            if (len(path_parts) != 2): break
            # process index file
            process_index_file(path_parts[0], filename, genome)
            
            # process collections file
            #collection_file = os.path.abspath(os.path.join(path_parts[0], 'collection.txt'))
            #cpath = Path(collection_file)
            #if (cpath.exists() and cpath.is_file()):
                # we also have a collection file      
            #    process_collection_file(path_parts[0], cpath)
        except Exception as e:
            print("Error accessing file: {0}".format(e))

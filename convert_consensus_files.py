import os
import argparse
import sys
from io_utility import *


# -----------------------------------------------------------------------------
# Reformat Files
# -----------------------------------------------------------------------------


def reformat_files(consensus_dir, output_dir):
    consensus_files = get_file_list(consensus_dir, ".fasta")
    for f in consensus_files:
        segment = os.path.basename(f).split(".")[1]
        with open(f, 'rb') as infile:
            for line in infile:
                if line[0] == ">":
                    line = line.strip('\n')
                    sample = line.split(" ")[0][1:]
                else:
                    line = line.strip('\n')
                    with open(output_dir+"/"+sample+"_consensus.fasta", 'a') as outfile:
                        outfile.write(">" + segment + "\n")
                        outfile.write(line + "\n")


def main():
    parser = argparse.ArgumentParser(description=("converts consensus files so that each file is a sample with all segment "
                                                  "consensus sequences"))

    parser.add_argument('consensus_dir', type=str, help='path to directory with the original consensus files')
    parser.add_argument('output_dir', type=str, help='path to the directory for the reformatted consensus files')

    if len(sys.argv) <= 1 or len(sys.argv) > 3:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    setup_output_dir(args.output_dir)
    reformat_files(args.consensus_dir, args.output_dir)


if __name__ == "__main__":
    main()

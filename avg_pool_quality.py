import pysam
import os
import fnmatch
import argparse
import sys
import csv
import numpy
from io_utility import *

# ---------------------------------------------------------------------------------------------
# Quality Assessment
# ---------------------------------------------------------------------------------------------


def get_quality(pacbio_dir, ref_seq, pools, output_dir):
    pacbio_bams = get_file_list(pacbio_dir, ".bam")
    segments = get_ref_info(ref_seq)
    for pool in pools:
        for segment in segments:
            quality_dict = dict()
            for bam in pacbio_bams:
                sample_id = "_".join(os.path.basename(bam).split("_")[0:2])
                if sample_id in pools[pool]:
                    bf = pysam.AlignmentFile(bam, "rb")
                    segment_length = len(segments[segment])
                    reads = bf.fetch(segment, 0, segment_length)
                    # holds the Read Objects
                    for read in reads:
                        for position in xrange(0, len(read.get_reference_positions())):
                            if position in quality_dict:
                                quality_dict[position].append(read.query_alignment_qualities[position])
                            else:
                                quality_dict[position] = [read.query_alignment_qualities[position]]
            write_output(pool, segment, quality_dict, output_dir)


def write_output(name, segment, quality_dict, output_dir):
    if not os.path.isfile(output_dir + name + "_avg_qual.csv"):
        with open(output_dir + name + "_avg_qual.csv", "a") as outfile:
            writer = csv.writer(outfile)
            writer.writerow(["pool", "segment", "position", "quality", "std", "num_reads"])
    with open(output_dir + name + "_avg_qual.csv", "a") as outfile:
            writer = csv.writer(outfile)
            for position in quality_dict:
                avg = numpy.mean(quality_dict[position])
                std = numpy.std(quality_dict[position])
                n = len(quality_dict[position])
                writer.writerow([name, segment, position, avg, std, n])


def main():
    parser = argparse.ArgumentParser(description=("asdasdasd"))

    parser.add_argument('pacbio_dir', type=str, help='path to the directory with the pacbio BAM files')
    parser.add_argument('ref_seq', type=str, help='path to reference sequence used in mapping')
    parser.add_argument('meta_data', type=str, help='path to the file with the pool info'),
    parser.add_argument('output_dir', type=str, help='path to the directory for the output files')

    if len(sys.argv) <= 3 or len(sys.argv) > 5:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    setup_output_dir(args.output_dir)
    pools = get_pool_info(args.meta_data)
    get_quality(args.pacbio_dir, args.ref_seq, pools, args.output_dir)


if __name__ == "__main__":
    main()

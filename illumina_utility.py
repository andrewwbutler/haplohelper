### find illumina reads that span all the positions of interest
from process_reads import *
import pysam


def illumina_haplotype(sample, variant_positions, segment, min_quality, output_dir):
    bf = pysam.AlignmentFile(sample.illumina_location, "rb")
    reads = bf.fetch(segment, 0, len(segment))
    filtered_reads = []
    for read in reads:
        # check that read passes filters
        read_obj = get_read_info(sample, read, segment)
        if illumina_filter(read_obj, min_quality, variant_positions):
            filtered_reads.append(read_obj)


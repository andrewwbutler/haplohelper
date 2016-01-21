from process_reads import *
import pysam


def illumina_haplotype(sample, variant_positions, segment, segment_length, min_quality, output_dir):
    bf = pysam.AlignmentFile(sample.illumina_location, "rb")
    reads = bf.fetch(segment)
    filtered_reads = []
    read_lengths = []
    for read in reads:
        read_lengths.append(len(read.query_alignment_sequence))
        # check that read passes filters
        read_obj = get_read_info(sample, read, segment)
        if illumina_filter(read_obj, min_quality, variant_positions):
            filtered_reads.append(read_obj)
    coverage_dict = dict()

    for read in filtered_reads:
        sequence = ""
        for position in variant_positions:
            sequence += read.sequence[int(position)-1]
        # ignore reads with insertions/deletions in positions of interest
        if "-" in sequence or len([c for c in sequence if c.islower()]) > 0:
            continue
        if sequence in coverage_dict:
            coverage_dict[sequence] += 1
        else:
            coverage_dict[sequence] = 1

    write_illumina_haplotypes(sample, segment, variant_positions, coverage_dict, output_dir)


def write_illumina_haplotypes(sample, segment, variant_positions, coverage_dict, output_dir):
    with open(output_dir + "/" + segment + "_illumina_linkage_info.csv", "a") as outfile:
        writer = csv.writer(outfile, delimiter=",")
        sample_name = sample.id.split("_")[0]
        if "Stock" in sample_name:
            day = "NA"
        else:
            day = sample.id.split("_")[1]

        writer.writerow(["sample", "day", "segment", "haplotype", "count", "freq", "variant_position", "nt"])

        for sequence, count in coverage_dict.items():
            idxs = [i for i, ltr in enumerate(sequence) if ltr != "N"]
            nts = [sequence[i] for i, ltr in enumerate(sequence) if ltr != "N"]
            for sequence2, count2 in coverage_dict():



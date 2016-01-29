from process_reads import *
import pysam


def illumina_haplotype(sample, variant_positions, segment, segment_length, min_quality, pb_haplotypes, output_dir):
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

    write_illumina_haplotypes(sample, segment, variant_positions, coverage_dict, pb_haplotypes, output_dir)


def write_illumina_haplotypes(sample, segment, variant_positions, coverage_dict, pb_haplotypes, output_dir):
    with open(output_dir + "/" + segment + "_illumina_linkage_info.csv", "a") as outfile:
        writer = csv.writer(outfile, delimiter=",")
        sample_name = sample.id.split("_")[0]
        if "Stock" in sample_name:
            day = "NA"
        else:
            day = sample.id.split("_")[1]

        il_haplotype_dict = dict()
        for sequence, count in coverage_dict.items():
            idxs = [i for i, ltr in enumerate(sequence) if ltr != "N"]
            nts = [sequence[i] for i, ltr in enumerate(sequence) if ltr != "N"]
            total_count = 0
            matching_count = 0
            # find all items in dictionary that have non N entries at idx specified, add to count
            for sequence2, count2 in coverage_dict.items():
                idxs2 = [i for i, ltr in enumerate(sequence2) if ltr != "N"]
                nts2 = [sequence2[i] for i, ltr in enumerate(sequence2) if ltr != "N"]
                if set(idxs) < set(idxs2) or idxs == idxs2:
                    total_count += count2
                    if nts == nts2:
                        matching_count += count2
            il_haplotype_dict[sequence] = (matching_count, total_count)

        pb_haplotype_dict = dict()
        for pb_haplotype, pb_count in pb_haplotypes:
            hap = ""
            for position in variant_positions:
                hap += pb_haplotype[int(position) - 1]
            pb_haplotype_dict[hap] = pb_count

        for pb_haplotype, pb_count in pb_haplotype_dict.items():
            i = 0
            for position in variant_positions:
                writer.writerow([sample_name, day, segment, pb_haplotype, pb_count, pb_haplotype, pb_count, 0, position, pb_haplotype[i]])
                i += 1
            for il_haplotype, il_count in il_haplotype_dict.items():
                idxs = [i for i, ltr in enumerate(il_haplotype) if ltr != "N"]
                if len(idxs) < 2:
                    continue
                keep = True
                for idx in idxs:
                    if il_haplotype[idx] != pb_haplotype[idx]:
                        keep = False
                        break
                if keep:
                    for idx in idxs:
                        writer.writerow([sample_name, day, segment, pb_haplotype, pb_count, il_haplotype, il_count[0], float(il_count[0])/il_count[1], variant_positions[idx], il_haplotype[idx]])

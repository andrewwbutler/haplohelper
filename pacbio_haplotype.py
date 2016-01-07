import pysam
import os
import argparse
import operator
from io_utility import *
from process_reads import *
import config
import csv

# ---------------------------------------------------------------------------------------------
# Haplotyping
# ---------------------------------------------------------------------------------------------


def gather_data(pacbio_dir, illumina_dir, consensus_dir, bino_filter, meta_data, samples):
    '''
    This function will create a list of Sample objects from the list of sample IDs provided.
    See process_reads.py for description of Sample objects.
    '''
    pacbio_bams = get_file_list(pacbio_dir, ".bam")
    illumina_snps = get_file_list(illumina_dir, "snplist.csv")
    consensus_files = get_file_list(consensus_dir, ".fasta")

    # create list of Sample objects from list of IDs given
    sample_list = []
    for bam in pacbio_bams:
        sample_id = "_".join(os.path.basename(bam).split("_")[0:2])
        if sample_id not in samples:
            continue
        snps = get_snps(illumina_snps, sample_id, bino_filter)
        if snps is None:
                print sample.id + " does not have an associated SNP list, skipping sample"
                continue
        consensus = get_consensus(consensus_files, sample_id)
        pools = get_pool_info(meta_data)
        for pool, sample in pools.iteritems():
            if sample_id in sample:
                sample_pool = pool
        sample_list.append(Sample(sample_id, bam, consensus, snps, sample_pool))
    return sample_list


def find_haplotypes(samples, segments, quality_dir, output_dir):
    for segment in segments:
        variant_positions = []
        all_haplotypes = dict()
        for sample in samples:
            bf = pysam.AlignmentFile(sample.read_location, "rb")
            for seg, sequence in sample.consensus.items():
                if seg == segment:
                    snps = []
                    for snp in sample.snps:
                        if snp[0] == segment:
                            snps.append(int(snp[1]))
                    if segment in sample.read_count_dict:
                        sample.read_count_dict[segment] += 1
                    else:
                        sample.read_count_dict[segment] = 1
                    segment_length = len(sequence)
                    reads = bf.fetch(segment, 0, segment_length)

                    for read in reads:
                        # check that read covers most of the segment
                        # check against illumina snp data
                        # check quality against pool quality
                        read_obj = get_read_info(sample, read, segment)
                        if skip_read(read_obj, segment, segment_length, snps, quality_dir, sample):
                            continue
                        sample.reads.append(read_obj)

                    haplotypes = dict()
                    for read in sample.reads:
                        haplotype = sequence
                        if len(read.mutations) == 0 and sequence not in haplotypes:
                            haplotypes[sequence] = 1
                        if sequence in haplotypes:
                            haplotypes[sequence] += 1
                        for mutation in read.mutations.items():
                            position = int(mutation[0]) - 1
                            haplotype = haplotype[0:position] + mutation[1] + haplotype[position+1:]
                            if int(mutation[0]) not in variant_positions:
                                variant_positions.append(int(mutation[0]))
                        if haplotype not in haplotypes:
                            haplotypes[haplotype] = 1
                        else:
                            haplotypes[haplotype] += 1

                    total = 0
                    for haplotype, count in haplotypes.items():
                        total += count

                    sorted_haplotypes = sorted(haplotypes.items(), key=operator.itemgetter(1), reverse=True)
                    # write_sample_haplotypes(sorted_haplotypes, output_dir)
                    all_haplotypes[sample.id] = sorted_haplotypes

        variant_positions = include_consensus_changes(samples, segment, variant_positions)
        variant_positions = sorted(variant_positions)
        # recheck reads so that they cover all variant positions

        for sample in samples:
            segment_length = len(sample.consensus[segment])
            for read in sample.reads:
                if skip_read(read, segment, segment_length, variant_positions, quality_dir, sample):
                    sample.reads.remove(read)

        for sample in samples:
            if len(sample.reads) > 0:
                write_haplotypes(sample, segment, all_haplotypes[sample.id], variant_positions, output_dir)
            else:
                print sample.id + " doesn't contain any reads that pass the filters"

        quit()


def include_consensus_changes(samples, segment, variant_positions):
    for sample1 in samples:
        s1 = sample1.consensus[segment]
        for sample2 in samples:
            s2 = sample2.consensus[segment]
            changes = [i for i in xrange(len(s1)) if s1[i] != s2[i]]
            for change in changes:
                if change not in variant_positions:
                    variant_positions.append(change)
    return variant_positions


def write_haplotypes(sample, segment, haplotypes, variant_positions, output_dir):
    with open(output_dir + "/haplotype_comp.csv", "a") as outfile:
        outfile.write(sample.id + "\n")
        writer = csv.writer(outfile, delimiter=",")
        consensus_positions = []
        for position in variant_positions:
            consensus_positions.append(sample.consensus[segment][int(position)-1])
        variant_positions.insert(0, "variant position")
        writer.writerow(variant_positions)
        consensus_positions.insert(0, "consensus")
        writer.writerow(consensus_positions)
        del variant_positions[0]

        total = 0
        for haplotype, count in haplotypes:
            total += count

        idx = 1
        for haplotype, count in haplotypes:
            row = ["haplotype_"+str(idx)]
            idx += 1
            for position in variant_positions:
                row.append(haplotype[int(position)-1])
            row.append(count)
            row.append(float(count)/total)
            writer.writerow(row)


def write_sample_haplotypes(sample, output_dir):
    with open(output_dir + "/haplotypes.txt", 'a') as outfile:
        x = 0
        for haplotype in sorted_haplotypes:
            x += 1
            outfile.write(">" + sample.id + "_" + str(x) + " " + str(haplotype[1]) + " " + str(float(haplotype[1])/total) + "\n")
            outfile.write(haplotype[0] + "\n")


def main():
    parser = argparse.ArgumentParser(description=("Python script for pacbio haplotyping. Argument parsing has been deprecated to the\
                                                configuration file. Help descriptions remain accurate"))

    parser.add_argument('pacbio_dir', type=str, help='path to the directory with the pacbio BAM files')
    parser.add_argument('illumina_dir', type=str, help='path to the directory with the illumina snp data')
    parser.add_argument('consensus_dir', type=str, help='path to per sample consensus data')
    parser.add_argument('quality_dir', type=str, help='path to the directory with pool quality info')
    parser.add_argument('barcode_meta_data', type=str, help='path to the file with the pool info')

    parser.add_argument('output_dir', type=str, help='path to the directory for the output files')
    parser.add_argument("-b", "--bino_filter", dest="bino_filter", action="store_false",
                        help='whether illumina snps must pass the distribution check, default is true')

    setup_output_dir(config.output_dir)
    sample_ids = get_sample_ids(config.samples_to_compare, config.sample_meta_data)
    samples = gather_data(config.pacbio_dir, config.illumina_dir, config.consensus_dir, config.bino_filter, config.barcode_meta_data,
                          sample_ids)
    segments = get_segments(samples, config.segments)
    find_haplotypes(samples, segments, config.quality_dir, config.output_dir)

if __name__ == "__main__":
    main()

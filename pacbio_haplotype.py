import pysam
import os
import argparse
import operator
from io_utility import *
from process_reads import *
from illumina_utility import *
import config
import csv

# ---------------------------------------------------------------------------------------------
# Haplotyping
# ---------------------------------------------------------------------------------------------


def gather_data(pacbio_dir, illumina_dir, illumina_bam_dir, consensus_dir, bino_filter, meta_data, samples, min_snp_freq, exclude_positions):
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
        sample_id = sample_id.split(".")[0]
        if sample_id.split("_")[1] == "A":
            sample_id = sample_id.split("_")[0]
        if sample_id not in samples:
            continue
        snps = get_snps(illumina_snps, sample_id, bino_filter, min_snp_freq, exclude_positions)
        if snps is None:
            print sample_id + " does not have an associated SNP list, skipping sample"
            continue
        consensus = get_consensus(consensus_files, sample_id)
        pools = get_pool_info(meta_data)
        illumina_location = get_illumina_location(illumina_bam_dir, sample_id)
        for pool, sample in pools.iteritems():
            if sample_id in sample:
                sample_pool = pool
        sample_list.append(Sample(sample_id, bam, illumina_location, consensus, snps, sample_pool))
    return sample_list


def find_haplotypes(samples, segments, quality_dir, illumina_dir, min_quality, output_dir):
    '''
    Main function for determining haplotypes
    '''
    for segment in segments:
        variant_positions = []
        all_haplotypes = dict()
        for sample in samples:
            print "Processing " + sample.id
            sequence = sample.consensus[segment]
            snps = get_segment_snps(sample, segment)
            segment_length = len(sequence)
            bf = pysam.AlignmentFile(sample.read_location, "rb")
            reads = bf.fetch(segment, 0, segment_length)

            for read in reads:
                # check that read passes filters, store those that do in sample.reads
                read_obj = get_read_info(sample, read, segment)
                if skip_read(read_obj, segment, segment_length, snps, quality_dir, sample):
                    continue
                sample.reads.append(read_obj)

            haplotypes = dict()
            for read in sample.reads:
                haplotype = sequence
                # corresponds to the consensus sequence
                if len(read.mutations) == 0 and sequence not in haplotypes:
                    haplotypes[sequence] = 1
                if sequence in haplotypes:
                    haplotypes[sequence] += 1
                for mutation in read.mutations.items():
                    if mutation[0] in snps:
                        position = int(mutation[0]) - 1
                        haplotype = haplotype[0:position] + mutation[1] + haplotype[position+1:]
                        if int(mutation[0]) not in variant_positions:
                            variant_positions.append(int(mutation[0]))
                if haplotype not in haplotypes:
                    haplotypes[haplotype] = 1
                else:
                    haplotypes[haplotype] += 1

            sorted_haplotypes = sorted(haplotypes.items(), key=operator.itemgetter(1), reverse=True)
            # write_sample_haplotypes(sorted_haplotypes, output_dir)
            all_haplotypes[sample.id] = sorted_haplotypes

        variant_positions = include_consensus_changes(samples, segment, variant_positions)
        variant_positions = sorted(variant_positions)

        with open(output_dir + "/" + segment + "_illumina_linkage_info.csv", "w") as outfile:
            writer = csv.writer(outfile, delimiter=",")
            writer.writerow(["sample", "day", "segment", "pb_haplotype", "pb_count", "il_haplotype", "il_count", "freq", "variant_position", "nt"])

        # recheck reads so that they cover all variant positions
        for sample in samples:
            if "Stock" not in sample.id:
                illumina_haplotype(sample, variant_positions, segment, segment_length, min_quality, all_haplotypes[sample.id], output_dir)
            segment_length = len(sample.consensus[segment])
            for read in sample.reads:
                if skip_read(read, segment, segment_length, variant_positions, quality_dir, sample):
                    sample.reads.remove(read)

        with open(output_dir + "/" + segment + "_tidy_haplotypes.csv", "w") as outfile:
            writer = csv.writer(outfile, delimiter=",")
            writer.writerow(["sample", "day", "segment", "haplotype", "count", "freq"])

        with open(output_dir + "/" + segment + "_tidier_haplotypes.csv", "w") as outfile:
            writer = csv.writer(outfile, delimiter=",")
            writer.writerow(["sample", "day", "segment", "haplotype", "count", "freq", "variant_position", "nt"])

        with open(output_dir + "/" + segment + "_illumina_minor_variants.csv", "w") as outfile:
            writer = csv.writer(outfile, delimiter=",")
            writer.writerow(["sample", "day", "segment", "position", "nt", "freq"])

        sample_names = []
        for sample in samples:
            if sample.id.split("_")[0] not in sample_names:
                sample_names.append(sample.id.split("_")[0])
            if len(sample.reads) > 0:
                write_segment_haplotypes(sample, segment, all_haplotypes[sample.id], variant_positions, illumina_dir, output_dir)
            else:
                print sample.id + " doesn't contain any reads that pass the filters"

        rank_haplotypes(segment, output_dir, variant_positions)
        check_for_takeover(output_dir, segment, samples, sample_names, all_haplotypes, variant_positions, illumina_dir)


def count_haplotypes(haplotypes):
    '''
    Given list of tuples with the haplotypes as first item and the number of reads corresponding
    to that haplotype as second item, calculate the total number of reads
    '''
    total = 0
    for haplotype, count in haplotypes:
        total += int(count)
    return total


def include_consensus_changes(samples, segment, variant_positions):
    '''
    Includes positions that are different between the consensus sequences of the samples
    in the variant_position list
    '''
    for sample1 in samples:
        s1 = sample1.consensus[segment]
        snps1 = get_segment_snps(sample1, segment)
        for sample2 in samples:
            snps2 = get_segment_snps(sample2, segment)
            s2 = sample2.consensus[segment]
            changes = [i for i in xrange(len(s1)) if s1[i] != s2[i]]
            for change in changes:
                if change not in variant_positions and (change in snps1 or change in snps2):
                    variant_positions.append(change)
    return variant_positions


def rank_haplotypes(segment, output_dir, variant_positions):
    '''
    Write a ranked list of all the haplotypes present in all the samples provided
    '''
    with open(output_dir + "/" + segment + "_haplotype_comp.csv", 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        haplotypes = dict()
        for row in reader:
            if len(row) != 0 and row[0].find("haplotype") != -1:
                haplotype = ""
                for i in xrange(1, len(variant_positions)+1):
                    haplotype += row[i]
                if haplotype in haplotypes:
                    haplotypes[haplotype] += int(row[len(variant_positions) + 1])
                else:
                    haplotypes[haplotype] = int(row[len(variant_positions) + 1])

        sorted_haplotypes = sorted(haplotypes.items(), key=operator.itemgetter(1), reverse=True)

    with open(output_dir + "/" + segment + "_ranked_haplotypes.csv", "w") as outfile:
        writer = csv.writer(outfile, delimiter=",")
        variant_positions.insert(0, "variant position")
        variant_positions.append("count")
        variant_positions.append("frequency")
        #writer.writerow(variant_positions)
        writer.writerow(["rank", "haplotype", "count", "frequency"])
        del variant_positions[0]
        del variant_positions[len(variant_positions) - 1]
        del variant_positions[len(variant_positions) - 1]

        total = count_haplotypes(sorted_haplotypes)
        idx = 1
        for haplotype in sorted_haplotypes:
            seq = [str(idx)]
            idx += 1
            hap_seq = ""
            for i in xrange(0, len(haplotype[0])):
                hap_seq += haplotype[0][i]
            seq.append(hap_seq)
            seq.append(haplotype[1])
            seq.append(str(float(haplotype[1])/total))
            seq.append(variant_positions)
            writer.writerow(seq)


def check_for_takeover(output_dir, segment, samples, sample_names, all_haplotypes, variant_positions, illumina_dir):
    '''
    Within a given sample (multiple days), check if there is a change in the major haplotype and whether the positions that change
    are confirmed in the illumina data.
    '''
    for name in sample_names:
        major_hap = ""
        for sample in samples:
            if not all_haplotypes[sample.id]:
                continue
            if sample.id.split("_")[0] == name:
                hap = ""
                for position in variant_positions:
                    hap += all_haplotypes[sample.id][0][0][int(position) - 1]
                illumina_hap = "".join(get_illumina_snv(illumina_dir, sample.id, segment, variant_positions)[2])
                if major_hap == "":
                    major_hap = hap
                elif major_hap != hap:
                    if hap == illumina_hap:
                        with open(output_dir + "/" + segment + "_haplotype_changes.txt", "a") as outfile:
                            outfile.write("Sample " + sample.id + " major haplotype changed from " + major_hap + " to " + hap +
                                          "---- confirmed by illumina" + "\n")

                    else:
                        with open(output_dir + "/" + segment + "_haplotype_changes.txt", "a") as outfile:
                            outfile.write("Sample " + sample.id + " major haplotype changed from " + major_hap + " to " + hap +
                                          "---- NOT confirmed by illumina" + "\n")
                    major_hap = hap
                else:
                    continue


def write_segment_haplotypes(sample, segment, haplotypes, variant_positions, illumina_dir, output_dir):
    '''
    Produces the main output file (one per segment) with the haplotype lists
    '''
    with open(output_dir + "/" + segment + "_haplotype_comp.csv", "a") as outfile:
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

        snvs = get_illumina_snv(illumina_dir, sample.id, segment, variant_positions)
        nts = snvs[0]
        nts.insert(0, "SNV")
        freqs = snvs[1]
        freqs.insert(0, "minor frequency")
        writer.writerow(nts)
        writer.writerow(freqs)

        total = count_haplotypes(haplotypes)

        idx = 1
        for haplotype, count in haplotypes:
            row = ["haplotype_"+str(idx)]
            idx += 1
            for position in variant_positions:
                row.append(haplotype[int(position)-1])
            row.append(count)
            row.append(float(count)/total)
            writer.writerow(row)
        outfile.write("\n")

    with open(output_dir + "/" + segment + "_tidy_haplotypes.csv", "a") as outfile:
        writer = csv.writer(outfile, delimiter=",")
        sample_name = sample.id.split("_")[0]
        if "Stock" in sample_name:
            day = "NA"
        else:
            day = sample.id.split("_")[1]
        total = count_haplotypes(haplotypes)
        for haplotype, count in haplotypes:
            sequence = ""
            for position in variant_positions:
                if sequence != "":
                    sequence += ", " + haplotype[int(position) - 1] + " - " + str(position)
                else:
                    sequence += haplotype[int(position) - 1] + " - " + str(position)
            freq = float(count)/total
            writer.writerow([sample_name, day, segment, sequence, count, freq])

    with open(output_dir + "/" + segment + "_tidier_haplotypes.csv", "a") as outfile:
        writer = csv.writer(outfile, delimiter=",")
        sample_name = sample.id.split("_")[0]
        if "Stock" in sample_name:
            day = "NA"
        else:
            day = sample.id.split("_")[1]
        total = count_haplotypes(haplotypes)
        idx = 0
        for haplotype, count in haplotypes:
            idx += 1
            hap = ""
            for position in variant_positions:
                hap += haplotype[int(position) - 1]
            for position in variant_positions:
                nt = haplotype[int(position) - 1]
                freq = float(count)/total
                writer.writerow([sample_name, day, segment, hap, count, freq, position, nt])

    with open(output_dir + "/" + segment + "_illumina_minor_variants.csv", "a") as outfile:
        writer = csv.writer(outfile, delimiter=",")
        if "Stock" in sample_name:
            day = "NA"
        else:
            day = sample.id.split("_")[1]
        sample_name = sample.id.split("_")[0]
        snvs = get_illumina_snv(illumina_dir, sample.id, segment, variant_positions)
        for i in xrange(0, len(snvs[0])):
            writer.writerow([sample_name, day, segment, variant_positions[i], snvs[0][i], snvs[1][i]])
            writer.writerow([sample_name, day, segment, variant_positions[i], snvs[2][i], snvs[3][i]])


def write_sample_haplotypes(sample, output_dir):
    '''
    Optional function to write all the full haplotypes within a single sample.
    '''
    with open(output_dir + "/haplotypes.txt", 'a') as outfile:
        x = 0
        for haplotype in sorted_haplotypes:
            x += 1
            outfile.write(">" + sample.id + "_" + str(x) + " " + str(haplotype[1]) + " " + str(float(haplotype[1])/total) + "\n")
            outfile.write(haplotype[0] + "\n")


def main():
    parser = argparse.ArgumentParser(description=("Python script for pacbio haplotyping. Argument parsing has been deprecated to the\
                                                configuration file. Please edit settings there. Help descriptions remain accurate"))

    parser.add_argument('pacbio_dir', type=str, help='path to the directory with the pacbio BAM files')
    parser.add_argument('illumina_dir', type=str, help='path to the directory with the illumina snp data')
    parser.add_argument('consensus_dir', type=str, help='path to per sample consensus data')
    parser.add_argument('quality_dir', type=str, help='path to the directory with pool quality info')
    parser.add_argument('barcode_meta_data', type=str, help='path to the file with the barcoding info')
    parser.add_argument('sample_meta_data', type=str, help='path to the file with the pool info')
    parser.add_argument('output_dir', type=str, help='path to the directory for the output files')
    parser.add_argument('segments', type=str, help="list of segments to analyze")
    parser.add_argument('samples_to_compare', type=str, help="list of samples to analyze or filters to select appropriate samples")

    parser.add_argument("-b", "--bino_filter", dest="bino_filter", action="store_false",
                        help='whether illumina snps must pass the distribution check, default is true')

    if len(sys.argv) > 1:
        parser.print_help()
        sys.exit(1)

    setup_output_dir(config.output_dir)
    sample_ids = get_sample_ids(config.samples_to_compare, config.sample_meta_data)
    samples = gather_data(config.pacbio_dir, config.illumina_dir, config.illumina_bam_dir, config.consensus_dir, config.bino_filter, config.barcode_meta_data, sample_ids, config.min_snp_freq, config.exclude_positions)
    segments = get_segments(samples, config.segments)
    find_haplotypes(samples, segments, config.quality_dir, config.illumina_dir, config.min_quality, config.output_dir)

if __name__ == "__main__":
    main()

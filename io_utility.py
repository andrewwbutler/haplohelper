import os
import sys
import fnmatch
import csv

# --------------------------------------------------------------------------------------------
# General Functions
# --------------------------------------------------------------------------------------------


def query_yes_no(question, default="yes"):
    '''
    Ask a yes/no question via raw_input() and decide whether to
        quit or continue.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).
    '''
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            if valid[choice]:
                return
            else:
                quit()
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")


def get_file_list(dir, suffix):
    '''
        Gathers all files in provided directory with the given extension and returns a list of paths to the files
    '''
    files = []
    for file in os.listdir(dir):
        if fnmatch.fnmatch(file, '*'+suffix):
            files.append(dir+file)
    return files


def setup_output_dir(output_dir):
    '''
        Creates the output directory, if it already exists it will prompt to delete all files within
    '''
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    else:
        #query_yes_no("Warning: Directory exists. All files within will be deleted. Continue?")
        fileList = os.listdir(output_dir)
        for fileName in fileList:
            os.remove(output_dir+"/"+fileName)

# --------------------------------------------------------------------------------------------
# Data Specific IO
# --------------------------------------------------------------------------------------------


def get_snps(snplists, sample, bino_filter, min_freq):
    '''
    Return a list of tuples where tuple[0] is the segment ID and tuple[1] is the snp position
    indexed from 1.
    '''
    snps = []
    for snplist in snplists:
        if snplist.find(sample) != -1:
            with open(snplist, 'rU') as csvfile:
                reader = csv.reader(csvfile, delimiter=',', quotechar='|')
                next(reader, None)
                for row in reader:
                    if float(row[7]) < min_freq:
                        continue
                    if bino_filter:
                        if row[3] == "TRUE" or row[3] == "True":
                            snp = (row[1], row[2])
                            snps.append(snp)
                    else:
                        snp = (row[1], row[2])
                        snps.append(snp)
            return snps


def get_consensus(consensus_files, sample):
    '''
    Return a dictionary to be stored in the Sample object where the key is the segment
    and the value is the corresponding consensus sequence.
    '''
    consensus_dict = dict()
    for consensus in consensus_files:
        if consensus.find(sample) != -1:
            with open(consensus, 'rb') as infile:
                for line in infile:
                    if line[0] == ">":
                        line = line.strip('\n')
                        segment = line.split(" ")[0][1:]
                        consensus_dict[segment] = ""
                    else:
                        line = line.strip('\n')
                        consensus_dict[segment] += line
            return consensus_dict


def get_ref_info(ref_seq):
    '''
    Return a list of all the segments present in the reference sequence.
    '''
    segments = dict()
    with open(ref_seq, "rU") as infile:
        current_segment = ""
        for line in infile:
            if line[0] == ">":
                line = line.strip('\n')
                segments[line[1:]] = ""
                current_segment = line[1:]
            else:
                line = line.strip('\n')
                segments[current_segment] += line
    return segments


def get_pool_info(meta_data):
    '''
    Given meta data csv with sample info as first entry and pool info as third entry,
    return a dictionary with pool as the key and list of corresponding samples as value.
    '''
    pool_info = dict()
    with open(meta_data, 'rU') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        next(reader, None)
        for row in reader:
            if row[2] in pool_info:
                pool_info[row[2]].append(row[0])
            else:
                pool_info[row[2]] = [row[0]]
    return pool_info


def build_quality_dict(quality_file, segment, snps):
    '''
    For a given quality file and segment, build a quality dictionary with key as position (index 1) and the
    value as the minimum accepted quality (mean-1sd).
    '''
    quality_dict = dict()
    with open(quality_file, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='l')
        next(reader, None)
        for row in reader:
            if row[1] == segment and (int(row[2])+1) in snps:
                quality_dict[int(row[2])+1] = float(row[3]) - float(row[4])
    return quality_dict


def get_segments(samples, segments):
    '''
    For a collection of samples return either a list of segments belonging to every sample
    or confirm that the list of segments passed are present in every sample.
    '''
    if segments[0] == "ALL":
        segments = []
        sample1 = samples[0]
        for segment, sequence in sample1.consensus.items():
            include = True
            for sample2 in samples:
                if segment not in sample2.consensus:
                    include = False
            if segment not in segments and include:
                segments.append(segment)
    else:
        for segment in segments:
            for sample in samples:
                contains = False
                for seg, seq in sample.consensus.items():
                    if seg in segments:
                        contains = True
                if not contains:
                    print "Sample " + sample.id + " does not contain " + segment + ". Removing from segment list."
                    segments.remove(segment)
    return segments


def get_sample_ids(samples_to_compare, meta_data):
    '''
    Enable user to easily filter which samples to include in the analysis. Can explicitly provide a list
    or can provide filters (day, generation, exposure status, exposure type) to generate the list of
    sample IDs that is returned.
    '''
    if samples_to_compare.samples and samples_to_compare.samples[0] != "ALL" and "_" in samples_to_compare.samples[0]:
        return samples_to_compare.samples
    else:
        samples = []
        sample_names = []
        if samples_to_compare.samples and samples_to_compare.samples[0] != "ALL":
            for sample in samples_to_compare.samples:
                sample_names.append(sample.split("_")[0])

        with open(meta_data, 'rU') as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='|')
            next(reader, None)
            for row in reader:
                if row[0] not in sample_names and sample_names:
                    continue
                if int(row[1]) not in samples_to_compare.days and samples_to_compare.days:
                    continue
                if row[5] not in samples_to_compare.generation and samples_to_compare.generation:
                    continue
                if row[9] not in samples_to_compare.prev_exposure_specific and samples_to_compare.prev_exposure_specific:
                    continue
                if row[12] not in samples_to_compare.prev_exposure_binary and samples_to_compare.prev_exposure_binary:
                    continue
                samples.append(row[0] + "_" + row[1])
    if samples_to_compare.samples and samples_to_compare.samples[0] == "ALL":
        samples.append("StockCal09")
    return samples


def get_illumina_snv(illumina_dir, sample, segment, variant_positions):
    '''
    Return the snv and frequency from the illumina data for each variant position
    '''
    snplists = get_file_list(illumina_dir, ".csv")
    snv = []
    frequency = []
    major = []
    major_freq = []
    for snplist in snplists:
        if snplist.find(sample) != -1:
            with open(snplist, 'rU') as csvfile:
                reader = csv.reader(csvfile, delimiter=',', quotechar='|')
                next(reader, None)
                for row in reader:
                    if int(row[2]) in variant_positions and row[1] == segment:
                        snv.append(row[6])
                        frequency.append(row[7])
                        major.append(row[4])
                        major_freq.append(row[5])
    return [snv, frequency, major, major_freq]


def get_transmission_pairs(meta_data, pb_meta_data):
    '''
    Returns a dictionary of all the f0/f1 pairs with f0 as the key, f1 as value
    '''
    samples = []
    with open(pb_meta_data, 'rU') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        next(reader, None)
        for row in reader:
            if row[0].split("_")[0] not in samples:
                samples.append(row[0].split("_")[0])

    pairs = dict()
    with open(meta_data, 'rU') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        next(reader, None)
        for row in reader:
            if row[5] == "F0" and row[0] in samples:
                pairs[row[0]] = row[7]
    return pairs

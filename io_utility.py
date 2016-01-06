import os
import sys
import fnmatch
import csv
import re

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
        query_yes_no("Warning: Directory exists. All files within will be deleted. Continue?")
        fileList = os.listdir(output_dir)
        for fileName in fileList:
            os.remove(output_dir+"/"+fileName)

# --------------------------------------------------------------------------------------------
# Data Specific IO
# --------------------------------------------------------------------------------------------


def get_snps(snplists, sample, bino_filter):
    snps = []
    for snplist in snplists:
        if snplist.find(sample) != -1:
            with open(snplist, 'rb') as csvfile:
                reader = csv.reader(csvfile, delimiter=',', quotechar='|')
                next(reader, None)
                for row in reader:
                    if bino_filter:
                        if row[3] == "True":
                            snp = (row[1], row[2])
                            snps.append(snp)
                    else:
                        snp = (row[1], row[2])
                        snps.append(snp)
            return snps


def get_consensus(consensus_files, sample):
    consensus_dict = dict()
    for consensus in consensus_files:
        if consensus.find(sample) != -1:
            with open(consensus, 'rb') as infile:
                for line in infile:
                    if line[0] == ">":
                        line = line.strip('\n')
                        sample = line.split(" ")[0][1:]
                        consensus_dict[sample] = ""
                    else:
                        line = line.strip('\n')
                        consensus_dict[sample] += line
            return consensus_dict


def get_ref_info(ref_seq):
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
    return(segments)


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


def build_quality_dict(quality_file, segment):
    '''
    For a given quality file and segment, build a quality dictionary with key as position and the
    value as the minimum accepted quality (mean-1sd)
    '''
    quality_dict = dict()
    with open(quality_file, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='l')
        next(reader, None)
        for row in reader:
            if row[1] == segment:
                quality_dict[row[2]] = float(row[3]) - float(row[4])
    return quality_dict


def get_sample_ids(samples_to_compare, meta_data):
    if samples_to_compare.samples:
        return samples_to_compare.samples
    else:
        samples = []
        with open(meta_data, 'rU') as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='|')
            next(reader, None)
            for row in reader:
                if int(row[1]) not in samples_to_compare.days and samples_to_compare.days:
                    continue
                if row[5] not in samples_to_compare.generation and samples_to_compare.generation:
                    continue
                if row[9] not in samples_to_compare.prev_exposure_specific and samples_to_compare.prev_exposure_specific:
                    continue
                if row[12] not in samples_to_compare.prev_exposure_binary and samples_to_compare.prev_exposure_binary:
                    continue
                samples.append(row[0] + "_" + row[1])
    return samples

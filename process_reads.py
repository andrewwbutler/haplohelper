from io_utility import *

# ---------------------------------------------------------------------------------------------
# Object Definitions
# ---------------------------------------------------------------------------------------------


class Sample(object):
    def __init__(self, sample_id, read_location, consensus, snps, pool):
        self.id = sample_id
        self.read_location = read_location
        self.consensus = consensus
        self.snps = snps
        self.reads = []
        self.read_count_dict = dict()
        self.pool = pool
        self.num_low_quality_reads = 0
        self.num_short_reads = 0

    def add_read(self, read):
        self.reads.append(read)


class Read(object):
    def __init__(self):
        self.sequence = ""
        self.mutations = dict()
        self.start_position = 0
        self.end_position = 0
        self.query_qualities = []

    def add_mutation(self, position, nt):
        self.mutations[position] = nt

    def add_nt(self, nt):
        self.sequence += nt


def skip_read(read, segment, segment_length, snps, quality_dir, sample):
    '''
    Determine whether or not to include a read in further analysis. Must pass
    a quality check and span from the first snp of interest to the last
    '''
    quality_files = get_file_list(quality_dir, ".csv")
    quality_dict = dict()
    for f in quality_files:
        if os.path.basename(f)[0] == sample.pool:
            quality_dict = build_quality_dict(f, segment, snps)
            break
    if not read_contains_snps(read, segment, snps, quality_dict, sample):
        return True
    if not pass_quality_check(read, snps, quality_dict, sample):
        return True
    return False


def pass_quality_check(read, snplist, quality_dict, sample):
    '''
    Performs the quality check - determines that the read is of sufficient
    quality at each required snp position
    '''
    start = read.start_position
    for snp in snplist:
        if read.query_qualities[int(snp) - 1 - start] < quality_dict[int(snp)]:
            sample.num_low_quality_reads += 1
            return False
    return True


def read_contains_snps(read, segment, snplist, quality_dict, sample):
    '''
    Performs the check to ensure that the read spans the first snp to the last
    '''
    start = read.start_position
    end = read.end_position
    for snp in snplist:
        if int(snp) < start + 1 or int(snp) > end + 1:
            sample.num_short_reads += 1
            return False
    return True


def get_read_info(sample, read, segment):
    '''
    Generates the read object for a given read. Also determines all substitution changes
    from the consensus.
    '''
    read_obj = Read()
    read_obj.start_position = read.get_reference_positions()[0]
    read_obj.end_position = read.get_reference_positions()[-1]
    read_obj.query_qualities = read.query_qualities
    sequence_idx = 0
    ref_idx = read.get_reference_positions()[0]
    ref_seq = sample.consensus[segment]

    for t in read.cigartuples:
        if t[0] == 0:
            for position in xrange(sequence_idx, sequence_idx + t[1]):
                if read.query_alignment_sequence[position] != ref_seq[ref_idx]:
                    read_obj.add_mutation(ref_idx+1, read.query_alignment_sequence[position])
                read_obj.add_nt(read.query_alignment_sequence[position])
                ref_idx += 1
                sequence_idx += 1
        if t[0] == 1:
            for position in xrange(sequence_idx, sequence_idx + t[1]):
                read_obj.add_nt(read.query_alignment_sequence[position].lower())
                sequence_idx += 1

        if t[0] == 2:
            for d in xrange(0, t[1]):
                read_obj.add_nt("-")
                ref_idx += 1

    # pad the end with Ns as necessary
    read_obj.add_nt("N"*(len(ref_seq) - ref_idx))
    return read_obj

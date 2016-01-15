# Argument list was getting a bit long for argparse. Use this configuration file
# to provide necessary arguments/settings

### File Locations

pacbio_dir = "/Users/abutler/Projects/flu_project/pacbio_sorted_bam/"
illumina_dir = "/Users/abutler/Projects/flu_project/illumina_snplists/"
consensus_dir = "/Users/abutler/Projects/flu_project/consensus/consensus_reformat/"
quality_dir = "/Users/abutler/Projects/flu_project/pacbio_quality/pool_quality/"
barcode_meta_data = "/Users/abutler/Projects/flu_project/pacbio_metadata_fixed.csv"
output_dir = "/Users/abutler/Projects/flu_project/pacbio_haplotypes/"
sample_meta_data = "/Users/abutler/Projects/flu_project/metadata.csv"


### Settings

# which segments to analyze for haplotypes, "ALL" will perform the evaluation for all segments common
# to every sample
segments = ["HA"]
# whether to include snps from the illumina data that don't pass the bino filter
# True -> don't include snps that don't pass
bino_filter = True
min_snp_freq = 0.01


### Samples of interest
class Sample_Specs(object):
    def __init__(self):
        self.samples = ["Gr5-3D", "Gr5-3C"]
        self.generation = []
        self.days = []
        self.prev_exposure_specific = []
        self.prev_exposure_binary = []


samples_to_compare = Sample_Specs()

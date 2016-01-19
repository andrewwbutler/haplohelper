### very hacky way to run all the transmission pairs and get plots for them.

from io_utility import *
import os
import argparse
import subprocess


def run_haplotyping(sample1, sample2):
    update_config_file(sample1, sample2)
    subprocess.call("python pacbio_haplotype.py", shell=True)
    subprocess.call("Rscript haplotype_visualization.R", shell=True)


def update_config_file(sample1, sample2):
    pattern = "self.samples"
    subst = "        self.samples = [\"" + sample1 + "\", \"" + sample2 + "\"]\n"
    new_file = open("new_config.py", 'w')
    old_file = open("config.py")
    for line in old_file:
        new_file.write(subst if pattern in line else line)

    new_file.close()
    old_file.close()
    os.remove("config.py")
    os.rename("new_config.py", "config.py")


def main():
    parser = argparse.ArgumentParser(description=(""))

    parser.add_argument('sample_meta_data', type=str, help='path to the file with the pool info')
    parser.add_argument('pb_meta_data', type=str, help='path to the file with the barcoding info')
    parser.add_argument('output_dir', type=str, help='path to the directory for the output files')

    if len(sys.argv) <= 2 or len(sys.argv) > 4:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    pairs = get_transmission_pairs(args.sample_meta_data, args.pb_meta_data)
    for f0, f1 in pairs.items():
        print "Comparing " + f0 + " vs " + f1
        run_haplotyping(f0, f1)
    #cleanup
    os.remove("Rplots.pdf")


if __name__ == "__main__":
    main()

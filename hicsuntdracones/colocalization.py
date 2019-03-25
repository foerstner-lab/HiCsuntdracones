from hicsuntdracones.colocalization_base import ColocalizationBase
import random
from collections import defaultdict
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt


class ColocalizationTester(ColocalizationBase):

    def __init__(self, matrix_file, gff_file, feature, bin_size,
                 output_prefix, number_of_subsamplings, margin,
                 flanks_only, remove_zeros):
        super().__init__(matrix_file, gff_file, feature, bin_size,
                 output_prefix, number_of_subsamplings, margin,
                 flanks_only, remove_zeros)
        self._interaction_matrix = None
        self._randomly_picked_bins = []
        self._feature_bin_countings = None
        self._random_bin_countings = None

    def generate_interaction_matrix(self):
        """Read the HiC matrix file and add each bin as further feature in the
                annotation.
                """
        print("- Reading HiC matrix file")
        self._interaction_matrix = self.add_matrix_bins_as_features\
            (self._matrix_file, self._features)

    def extract_random_bins(self):
        """For each feature of interest generate a equivalently long one from
        a random position and extract the ovelapping bins..
        """
        chrom_and_length = defaultdict(int)
        for genome_bin in self._interaction_matrix["Regions"]:
            chrom, length = genome_bin.split("-")
            if chrom_and_length[chrom] < int(length):
                chrom_and_length[chrom] = int(length)
        random.seed(1)
        for feature_of_interest in self._features.features_of_type(
            self._feature):
            # Perform the random selection x times:
            subsampling = self._number_of_subsamplings
            feature_length = (feature_of_interest.end -
                              feature_of_interest.start)
            while subsampling > 0:
                random_chrom = random.choice(list(chrom_and_length.keys()))
                if chrom_and_length[random_chrom] <= feature_length:
                    continue
                random_start = random.randint(
                    0, chrom_and_length[random_chrom] - feature_length)
                random_end = random_start + feature_length
                self._randomly_picked_bins.extend(
                    self._extract_feature_overlapping_bins(
                        random_start, random_end, random_chrom))
                subsampling -= 1

    def compile_interaction_valus_of_bins(self):
        self._feature_bin_countings = self._extract_selected_countings(
            self._feature_overlapping_bins, self._interaction_matrix)
        self._random_bin_countings = self._extract_selected_countings(
            self._randomly_picked_bins, self._interaction_matrix)
        if self._remove_zeros:
            self._feature_bin_countings = list(filter(
                lambda counting: counting != 0.0, self._feature_bin_countings))
            self._random_bin_countings = list(filter(
                lambda counting: counting != 0.0, self._random_bin_countings))

    def perform_t_test(self):
        with open(f"{self._output_prefix}_t-test_results.txt", "w") as output_fh:
            tstat, pvalue = stats.ttest_ind(
                self._feature_bin_countings, self._random_bin_countings,
                equal_var=False)
            output_fh.write(f"Margin: {self._margin}\n")
            output_fh.write(f"Feature of interest median: {np.median(self._feature_bin_countings)}\n")
            output_fh.write(f"Feature of interest mean: {np.mean(self._feature_bin_countings)}\n")
            output_fh.write(f"Feature of interest standard deviation: {np.std(self._feature_bin_countings)}\n")
            output_fh.write(f"Background median: {np.median(self._random_bin_countings)}\n")
            output_fh.write(f"Background mean: {np.mean(self._random_bin_countings)}\n")
            output_fh.write(f"Background standard deviation: {np.std(self._random_bin_countings)}\n")
            output_fh.write(f"Welch's t-test statistic: {tstat}\n")
            output_fh.write(f"Welch's t-test p-value: {pvalue}\n")

    def write_countings_for_file(self):
        with open(f"{self._output_prefix}_countings.txt", "w") as output_fh:
            feature_bin_countings_str = ", ".join([str(value) for value in self._feature_bin_countings])
            output_fh.write(f"Feature of interest ({self._feature}):\t{feature_bin_countings_str}\n")
            random_bin_countings_str = ", ".join([str(value) for value in self._random_bin_countings])
            output_fh.write(f"Background:\t{random_bin_countings_str}\n")

    def plot_distribution(self):
        plt.style.use('ggplot')
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.hist(np.array(self._feature_bin_countings), alpha=0.5, bins=200, label="feature", density=1)
        ax.hist(self._random_bin_countings, alpha=0.5, bins=200, label="background", density=1)
        fig.tight_layout()
        ax.legend(loc='upper center')
        plt.savefig(f"{self._output_prefix}_histograms.pdf")

def analyse_colocalization(args):
    colocalization_tester = ColocalizationTester(
        args.matrix_file, args.gff_file, args.feature, args.bin_size,
        args.output_prefix, args.number_of_subsamplings, args.margin,
        args.flanks_only, args.remove_zeros)
    colocalization_tester.read_gff_file()
    colocalization_tester.generate_interaction_matrix()
    colocalization_tester.extract_features_overlapping_bins()
    colocalization_tester.extract_random_bins()
    colocalization_tester.compile_interaction_valus_of_bins()
    colocalization_tester.perform_t_test()
    colocalization_tester.write_countings_for_file()
    colocalization_tester.plot_distribution()


def remove_position_information(name_with_pos_info):
    # Return just the chromosome part without the exact window
    # location
    return "-".join(name_with_pos_info.split("-")[:-1])

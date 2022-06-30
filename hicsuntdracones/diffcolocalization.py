from hicsuntdracones.colocalization_base import ColocalizationBase
import random
from collections import defaultdict
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt


class DiffColocalizationTester(ColocalizationBase):

    def __init__(self, matrix_file, other_matrix_file, gff_file, feature, bin_size,
                 output_prefix, number_of_subsamplings, margin,
                 flanks_only, remove_zeros):
        self._other_matrix_file = other_matrix_file
        super().__init__(matrix_file, gff_file, feature, bin_size,
                         output_prefix, number_of_subsamplings, margin,
                         flanks_only, remove_zeros)
        self._interaction_matrix_1 = None
        self._interaction_matrix_2 = None
        self._randomly_picked_bin_lists = []
        self._feature_bin_countings_1 = None
        self._feature_bin_countings_2 = None
        self._random_bin_countings_lists_1 = None
        self._random_bin_countings_lists_2 = None
        self._interaction_ratios_1 = None
        self._interaction_ratios_2 = None
        self._interaction_ratios_1_clean = None
        self._interaction_ratios_2_clean = None

    def generate_both_interaction_matrices(self):
        import copy
        """Read the HiC matrix file and add each bin as further feature in the
                annotation.
                """
        print("- Reading first HiC matrix file")
        self._interaction_matrix_1 = self.add_matrix_bins_as_features\
            (self._matrix_file, copy.copy(self._features))
        print("- Reading second HiC matrix file")
        self._interaction_matrix_2 = self.add_matrix_bins_as_features\
            (self._other_matrix_file, copy.copy(self._features))

    def extract_random_bins(self):
        """
        Need to be done only for one of the matrix feature lists as we will
        use the same random bin for both matrices.
        For each feature of interest generate a equivalently long,
        temporary feature starting from a random position and extract
        the bins ovelapping with it.
        For each subsampling we store the list of all bin ids as a
        list. This grouping by subsampling rounds is essential!
        """
        # Will become a list of lists

        chrom_and_length = defaultdict(int)
        for genome_bin in self._interaction_matrix_1["Regions"]:
            chrom, length = genome_bin.split("-")
            if chrom_and_length[chrom] < int(length):
                chrom_and_length[chrom] = int(length)
        random.seed(1)
        subsampling = self._number_of_subsamplings
        while subsampling > 0:
            subsampling -= 1
            bins_of_cur_subsampling = []
            for feature_of_interest in self._features.features_of_type(
                self._feature):
                feature_length = (feature_of_interest.end -
                                  feature_of_interest.start)
                random_chrom = random.choice(list(chrom_and_length.keys()))
                if chrom_and_length[random_chrom] <= feature_length:
                    continue
                random_start = random.randint(
                    0, chrom_and_length[random_chrom] - feature_length)
                random_end = random_start + feature_length
                bins_of_cur_subsampling.extend(
                    self._extract_feature_overlapping_bins(
                        random_start, random_end, random_chrom))
            self._randomly_picked_bin_lists.append(bins_of_cur_subsampling)

    def compile_interaction_valus_of_bins(self):
        # Here lists of bin ids are converted in lists of interaction values
        self._feature_bin_countings_1 = self._extract_selected_countings(
            self._feature_overlapping_bins, self._interaction_matrix_1)
        self._feature_bin_countings_2 = self._extract_selected_countings(
            self._feature_overlapping_bins, self._interaction_matrix_2)

        # Here list of list are translated: list of list of bin ids
        # are converted into list of list of interation values
        self._random_bin_countings_lists_1 = [
            self._extract_selected_countings(
                bin_id_list, self._interaction_matrix_1)
            for bin_id_list in self._randomly_picked_bin_lists]
        self._random_bin_countings_lists_2 = [
            self._extract_selected_countings(
                bin_id_list, self._interaction_matrix_2)
            for bin_id_list in self._randomly_picked_bin_lists]

    def calculate_interaction_ratios(self):

        """
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Attention: Bins containing values of zero are excluded from
        the colocalization analysis as they reflect e.g. repetetive
        reagions to which no read can be uniquely aligned. This can
        result to different samples sizes during the comparison.
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        """

        if not self._remove_zeros:
            self._interaction_ratios_1 = [
                np.median(self._feature_bin_countings_1)
                / np.median(random_bin_countings)
                for random_bin_countings in self._random_bin_countings_lists_1]
            self._interaction_ratios_2 = [
                np.median(self._feature_bin_countings_2)
                / np.median(random_bin_countings)
                for random_bin_countings in self._random_bin_countings_lists_2]
        else:
            self._interaction_ratios_1 = [
                np.median(self._feature_bin_countings_1[
                              self._feature_bin_countings_1.nonzero()])
                / np.median(random_bin_countings[
                                random_bin_countings.nonzero()])
                for random_bin_countings in self._random_bin_countings_lists_1]
            self._interaction_ratios_2 = [
                np.median(self._feature_bin_countings_2[
                              self._feature_bin_countings_2.nonzero()])
                / np.median(random_bin_countings[
                                random_bin_countings.nonzero()])
                for random_bin_countings in self._random_bin_countings_lists_2]

    def perform_t_test(self):
        # The ratio will be inf if the median of the feature bins is >
        # 0 and the median or the random bins is = 0. If both are 0
        # the ratio will be nan.
        self._interaction_ratios_1_clean = list(filter(
            lambda value: not (np.isnan(value) or np.isinf(value)),
            self._interaction_ratios_1))
        self._interaction_ratios_2_clean = list(filter(
            lambda value: not (np.isnan(value) or np.isinf(value)),
            self._interaction_ratios_2))

        with open(f"{self._output_prefix}_t-test_results.txt", "w") as output_fh:
            tstat, pvalue = stats.ttest_ind(
                self._interaction_ratios_1_clean,
                self._interaction_ratios_2_clean,
                equal_var=False)
            output_fh.write(f"Margin: {self._margin}\n")
            output_fh.write(f"Number of subsamplings: {self._number_of_subsamplings}\n")
            len_diff_1_str = str(len(self._interaction_ratios_1) -
                                 len(self._interaction_ratios_1_clean))
            output_fh.write(f"Number of inf/nan removed set 1: {len_diff_1_str}\n")
            len_diff_2_str = str(len(self._interaction_ratios_2) -
                                 len(self._interaction_ratios_2_clean))
            output_fh.write(f"Number of inf/nan removed set 2: {len_diff_2_str}\n")
            output_fh.write(f"Matrix 1 ratio median: {np.median(self._interaction_ratios_1_clean)}\n")
            output_fh.write(f"Matrix 1 ratio mean: {np.mean(self._interaction_ratios_1_clean)}\n")
            output_fh.write(f"Matrix 1 ratio deviation: {np.std(self._interaction_ratios_1_clean)}\n")
            output_fh.write(f"Matrix 2 ratio median: {np.median(self._interaction_ratios_2_clean)}\n")
            output_fh.write(f"Matrix 2 ratio mean: {np.mean(self._interaction_ratios_2_clean)}\n")
            output_fh.write(f"Matrix 2 ratio deviation: {np.std(self._interaction_ratios_2_clean)}\n")
            output_fh.write(f"Welch's t-test statistic: {tstat}\n")
            output_fh.write(f"Welch's t-test p-value: {pvalue}\n")

    def write_countings_for_file(self):
        with open(f"{self._output_prefix}_countings.txt", "w") as output_fh:
            interaction_ratios_1_str = ", ".join([str(value) for value in self._interaction_ratios_1])
            output_fh.write(f"Interaction ratios 1:\t{interaction_ratios_1_str}\n")
            interaction_ratios_2_str = ", ".join([str(value) for value in self._interaction_ratios_2])
            output_fh.write(f"Interaction ratios 2:\t{interaction_ratios_2_str}\n")
            interaction_ratios_1_clean_str = ", ".join([str(value) for value in self._interaction_ratios_1_clean])
            output_fh.write(f"Interaction ratios 1 cleaned:\t{interaction_ratios_1_clean_str}\n")
            _interaction_ratios_2_clean_str = ", ".join([str(value) for value in self._interaction_ratios_2_clean])
            output_fh.write(f"Interaction ratios 2 cleaned:\t{_interaction_ratios_2_clean_str}\n")

    def plot_distribution(self):
        plt.style.use('ggplot')
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.hist(self._interaction_ratios_1_clean, alpha=0.5, bins=200, label="Interaction ratios 1", density=1)
        ax.hist(self._interaction_ratios_2_clean, alpha=0.5, bins=200, label="Interaction ratios 2", density=1)
        fig.tight_layout()
        ax.legend(loc='upper center')
        plt.savefig(f"{self._output_prefix}_histograms.pdf")


def analyse_diff_colocalization(args):
    diff_colocalization_tester = DiffColocalizationTester(
        args.matrix_file_1, args.matrix_file_2, args.gff_file, args.feature,
        args.bin_size, args.output_prefix, args.number_of_subsamplings,
        args.margin, args.flanks_only, args.remove_zeros)
    diff_colocalization_tester.read_gff_file()
    diff_colocalization_tester.generate_both_interaction_matrices()
    diff_colocalization_tester.extract_features_overlapping_bins()
    diff_colocalization_tester.extract_random_bins()
    diff_colocalization_tester.compile_interaction_valus_of_bins()
    diff_colocalization_tester.calculate_interaction_ratios()
    diff_colocalization_tester.perform_t_test()
    diff_colocalization_tester.write_countings_for_file()
    diff_colocalization_tester.plot_distribution()


# if __name__ == "__main__":
#     main()

import argparse
from collections import defaultdict
import itertools
import random
import gffutils
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import hicsuntdracones.hicmatrix


def analyse_diff_colocalization(args):
    diff_colocalization_tester = DiffColocalizationTester(
        args.matrix_file_1, args.matrix_file_2, args.gff_file, args.feature,
        args.bin_size, args.output_prefix, args.number_of_subsamplings,
        args.margin, args.flanks_only, args.remove_zeros)
    diff_colocalization_tester.read_gff_file()
    diff_colocalization_tester.add_matrix_bins_as_features()
    diff_colocalization_tester.extract_features_overlapping_bins()
    diff_colocalization_tester.extract_random_bins()
    diff_colocalization_tester.compile_interaction_valus_of_bins()
    diff_colocalization_tester.calculate_interaction_ratios()
    diff_colocalization_tester.perform_t_test()
    diff_colocalization_tester.write_countings_for_file()
    diff_colocalization_tester.plot_distribution()


class DiffColocalizationTester(object):

    def __init__(self, matrix_file_1, matrix_file_2, gff_file, feature,
                 bin_size, output_prefix, number_of_subsamplings, margin,
                 flanks_only, remove_zeros):
        self._matrix_file_1 = matrix_file_1
        self._matrix_file_2 = matrix_file_2
        self._gff_file = gff_file
        self._feature = feature
        self._bin_size = bin_size
        self._output_prefix = output_prefix
        self._number_of_subsamplings = number_of_subsamplings
        self._margin = margin
        self._flanks_only = flanks_only
        self._remove_zeros = remove_zeros

    def read_gff_file(self):
        """Read the GFF file twice and store in seperated feature list as for
        each HiC matrix bins will be added as further features. This
        makes it possible to query bins overlapping with features for
        each matrix seperately.

        """
        print("- Reading GFF file")
        # Same GGF file is read twice. Later bins from differen
        # matrices will be added to them.
        self._features_1 = gffutils.create_db(self._gff_file, ":memory:")
        self._features_2 = gffutils.create_db(self._gff_file, ":memory:")
        
    def add_matrix_bins_as_features(self):
        """Read the HiC matrix file and add each bin as further feature to the
        annotation features.
        """
        print("- Reading HiC matrix files")
        self._interaction_matrix_1 = self._add_matrix_bins_as_features(
            self._matrix_file_1, self._features_1)
        self._interaction_matrix_2 = self._add_matrix_bins_as_features(
            self._matrix_file_2, self._features_2)

    def _add_matrix_bins_as_features(self, matrix_file, features):
        hic_matrix = hicsuntdracones.hicmatrix.HiCMatrix(matrix_file)
        hic_matrix.normalize_by_columns_sum()
        interaction_matrix = hic_matrix.hic_matrix_df
        interaction_matrix.set_index(
            interaction_matrix["Regions"], inplace=True)
        for genome_bin in interaction_matrix["Regions"]:
            chrom_part = "-".join(genome_bin.split("-")[:-1])
            # As the bin counting starts with 0 but the gff starts at 1
            # we have to add 1 here to the position.
            # pos = int(genome_bin.split("-")[-1]) + 1
            pos = int(genome_bin.split("-")[-1])
            pos_adjusted = pos + 1
            features.update([gffutils.Feature(
                seqid=chrom_part,
                source='-',
                featuretype='HiC_bin',
                start=pos_adjusted,
                end=pos_adjusted+self._bin_size,
                strand="+",
                attributes="ID={}".format(genome_bin))])
        return interaction_matrix

    def extract_features_overlapping_bins(self):
        print("- Extracting bins that overlap given features")
        self._feature_overlapping_bins = []
        for feature_of_interest in self._features_1.features_of_type(
                self._feature):
            cur_overlapping_bins = self._extract_feature_overlapping_bins(
                feature_of_interest.start, feature_of_interest.end,
                feature_of_interest.seqid)
            self._feature_overlapping_bins.extend(cur_overlapping_bins)

    def _extract_feature_overlapping_bins(self, feature_start, feature_end,
                                          feature_seqid):
        if not self._flanks_only:
            return self._extract_feature_overlapping_bins_feature_and_flanks(
                feature_start, feature_end, feature_seqid)
        else:
            return self._extract_feature_overlapping_bins_flanks_only(
                feature_start, feature_end, feature_seqid)

    def _extract_feature_overlapping_bins_feature_and_flanks(
            self, feature_start, feature_end, feature_seqid):
        """

          Flank           Feature             Flank
        |----------=========================----------|

                             |
                             |
                             v
 
           Region in which overlapping bins are used
        |---------------------------------------------|

        """
        # The column name is the SeqID and the
        # position after correcting for difference between the
        # 0- to 1-based difference.
        start = feature_start - (self._margin * self._bin_size)
        if start < 1:
            start = 1
        end = feature_end + (self._margin * self._bin_size)
        return self._extract_overlapping_bins(start, end, feature_seqid)

    def _extract_feature_overlapping_bins_flanks_only(
            self, feature_start, feature_end, feature_seqid):
        overlapping_bins = []
        """
          Flank           Feature             Flank
        |----------=========================----------|

                             |
                             |
                             v
 
           Region in which overlapping bins are used
        |----------|                       |----------|
        """
        # Flank upstream of feature
        #   *Flank*           Feature             Flank
        # |----------=========================----------|
        start = feature_start - (self._margin * self._bin_size)
        if start < 1:
            start = 1
        end = feature_start
        overlapping_bins.extend(
            self._extract_overlapping_bins(start, end, feature_seqid))

        # Flank downstream of feature
        #   Flank           Feature             *Flank*
        # |----------=========================----------|
        start = feature_end
        end = feature_end + (self._margin * self._bin_size)
        overlapping_bins.extend(
            self._extract_overlapping_bins(start, end, feature_seqid))
        return overlapping_bins

    def _extract_overlapping_bins(self, start, end, feature_seqid):
        return [overlapping_bin.attributes["ID"][0]
                for overlapping_bin in self._features_1.region(
                        seqid=feature_seqid,
                        start=start,
                        end=end,
                        featuretype="HiC_bin")]

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
        self._randomly_picked_bin_lists = []
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
            for feature_of_interest in self._features_1.features_of_type(
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

    def _extract_selected_countings(self, bin_list, interaction_matrix):
        """Here the values of the interaction between two bins are extracted
        from a given matrix.
        """
        combined_interaction_countings = []
        for interaction_bin in bin_list:
            # Only count the interaction with bins on other
            # chromosomes and cores only
            other_bins = bin_list.copy()
            other_bins.remove(interaction_bin)
            chrom = interaction_bin.split("-")[0]
            #  Interchromosomal only i.e. have to be located on other
            #  chromosomes.
            other_bins = list(filter(
                lambda interaction_bin:
                not interaction_bin.startswith(chrom),
                other_bins))
            interaction_countings = interaction_matrix.ix[
                [interaction_bin], other_bins]
            combined_interaction_countings.extend(
                list(itertools.chain(*interaction_countings.values.tolist())))
        combined_interaction_countings = np.array(
            combined_interaction_countings)
        return combined_interaction_countings

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
        
        with open("{}_t-test_results.txt".format(
                self._output_prefix), "w") as output_fh:
            tstat, pvalue = stats.ttest_ind(
                self._interaction_ratios_1_clean,
                self._interaction_ratios_2_clean,
                equal_var=False)
            output_fh.write("Margin: {}\n".format(self._margin))
            output_fh.write("Number of subsamplings: {}\n".format(
                self._number_of_subsamplings))
            output_fh.write("Number of inf/nan removed set 1: {}\n".format(
                len(self._interaction_ratios_1) -
                len(self._interaction_ratios_1_clean)))
            output_fh.write("Number of inf/nan removed set 2: {}\n".format(
                len(self._interaction_ratios_2) -
                len(self._interaction_ratios_2_clean)))
            output_fh.write("Matrix 1 ratio median: {}\n".format(
                np.median(self._interaction_ratios_1_clean)))
            output_fh.write("Matrix 1 ratio mean: {}\n".format(
                np.mean(self._interaction_ratios_1_clean)))
            output_fh.write("Matrix 1 ratio deviation: "
                            "{}\n".format(np.std(
                                self._interaction_ratios_1_clean)))
            output_fh.write("Matrix 2 ratio median: {}\n".format(
                np.median(self._interaction_ratios_2_clean)))
            output_fh.write("Matrix 2 ratio mean: {}\n".format(
                np.mean(self._interaction_ratios_2_clean)))
            output_fh.write("Matrix 2 ratio deviation: "
                            "{}\n".format(np.std(
                                self._interaction_ratios_2_clean)))
            output_fh.write("Welch's t-test statistic: {}\n".format(tstat))
            output_fh.write("Welch's t-test p-value: {}\n".format(pvalue))

    def write_countings_for_file(self):
        with open("{}_countings.txt".format(self._output_prefix),
                  "w") as output_fh:
            output_fh.write("Interaction ratios 1:\t{}\n".format(", ".join(
                    [str(value) for value in self._interaction_ratios_1])))
            output_fh.write("Interaction ratios 2:\t{}\n".format(", ".join(
                    [str(value) for value in self._interaction_ratios_2])))
            output_fh.write("Interaction ratios 1 cleaned:\t{}\n".format(
                ", ".join(
                    [str(value) for value
                     in self._interaction_ratios_1_clean])))
            output_fh.write("Interaction ratios 2 cleaned:\t{}\n".format(
                ", ".join(
                    [str(value) for value
                     in self._interaction_ratios_2_clean])))

    def plot_distribution(self):
        plt.style.use('ggplot')
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.hist(
            self._interaction_ratios_1_clean, alpha=0.5, bins=200,
            label="Interaction ratios 1", normed=1)
        ax.hist(
            self._interaction_ratios_2_clean, alpha=0.5, bins=200,
            label="Interaction ratios 2", normed=1)
        fig.tight_layout()
        ax.legend(loc='upper center')
        plt.savefig("{}_histograms.pdf".format(self._output_prefix))


if __name__ == "__main__":
    main()

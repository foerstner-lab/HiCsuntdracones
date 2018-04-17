from collections import defaultdict
import itertools
from scipy import stats
import numpy as np
import gffutils
import random
import hicsuntdracones.hicmatrix
import matplotlib
matplotlib.use("Agg")
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt


def analyse_colocalization(args):
    colocalization_tester = ColocalizationTester(
        args.matrix_file, args.gff_file, args.feature, args.bin_size,
        args.output_prefix, args.number_of_subsamplings, args.margin,
        args.flanks_only, args.remove_zeros)
    colocalization_tester.read_gff_file()
    colocalization_tester.read_matrix_file_and_add_as_features()
    colocalization_tester.extract_features_overlapping_bins()
    colocalization_tester.extract_random_bins()
    colocalization_tester.compile_interaction_valus_of_bins()
    colocalization_tester.perform_t_test()
    colocalization_tester.write_countings_for_file()
    colocalization_tester.plot_distribution()


class ColocalizationTester(object):
    
    def __init__(self, matrix_file, gff_file, feature, bin_size,
                 output_prefix, number_of_subsamplings, margin,
                 flanks_only, remove_zeros):
        self._matrix_file = matrix_file
        self._gff_file = gff_file
        self._feature = feature
        self._bin_size = bin_size
        self._output_prefix = output_prefix
        self._number_of_subsamplings = number_of_subsamplings
        self._margin = margin
        self._flanks_only = flanks_only
        self._remove_zeros = remove_zeros

    def read_gff_file(self):
        print("- Reading GFF file")
        self._features = gffutils.create_db(self._gff_file, ":memory:")
        
    def read_matrix_file_and_add_as_features(self):
        """Read the HiC matrix file and add each bin as further feature in the
        annotation.
        """
        print("- Reading HiC matrix file")
        hic_matrix = hicsuntdracones.hicmatrix.HiCMatrix(self._matrix_file)
        hic_matrix.normalize_by_columns_sum()
        self._interaction_matrix = hic_matrix.hic_matrix_df
        # self._interaction_matrix["Regions"] = self._interaction_matrix[
        #     "Regions"].apply(remove_position_information)
        self._interaction_matrix.set_index(
            self._interaction_matrix["Regions"], inplace=True)

        for genome_bin in self._interaction_matrix["Regions"]:
            chrom_part = "-".join(genome_bin.split("-")[:-1])
            # As the bin counting starts with 0 but the gff starts a 1
            # we have to add 1 here to the position.
            # pos = int(genome_bin.split("-")[-1]) + 1
            pos = int(genome_bin.split("-")[-1])
            pos_adjusted = pos + 1
            self._features.update([gffutils.Feature(
                seqid=chrom_part,
                source='-',
                featuretype='HiC_bin',
                start=pos_adjusted,
                end=pos_adjusted+self._bin_size,
                strand="+",
                attributes="ID={}".format(genome_bin))])

    def extract_features_overlapping_bins(self):
        print("- Extracting bins that overlap given features")
        self._feature_overlapping_bins = []
        for feature_of_interest in self._features.features_of_type(
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
                for overlapping_bin in self._features.region(
                        seqid=feature_seqid,
                        start=start,
                        end=end,
                        featuretype="HiC_bin")]
    
    def extract_random_bins(self):
        """For each feature of interest generate a equivalently long one from
        a random position and extract the ovelapping bins..

        """
        self._randomly_picked_bins = []
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
            self._feature_overlapping_bins)
        self._random_bin_countings = self._extract_selected_countings(
           self._randomly_picked_bins)
        if self._remove_zeros:
            self._feature_bin_countings = list(filter(
                lambda counting: counting != 0.0, self._feature_bin_countings))
            self._random_bin_countings = list(filter(
                lambda counting: counting != 0.0, self._random_bin_countings))

    def _extract_selected_countings(self, bin_list):
        combined_countings = []
        for interaction_bin in bin_list:
            # Only count the interaction with bins on other
            # chromosoms
            other_bins = bin_list.copy()
            other_bins.remove(interaction_bin)
            chrom = interaction_bin.split("-")[0]
            #  Interchromosomal only i.e. have to be located on other
            #  chromosomes.
            other_bins = list(filter(
                lambda interaction_bin:
                not interaction_bin.startswith(chrom),
                other_bins))
            countings = self._interaction_matrix.ix[
                [interaction_bin], other_bins]
            combined_countings.extend(
                list(itertools.chain(*countings.values.tolist())))
        return np.array(combined_countings)
        
    def perform_t_test(self):
        with open("{}_t-test_results.txt".format(
                self._output_prefix), "w") as output_fh:
            tstat, pvalue = stats.ttest_ind(
                self._feature_bin_countings, self._random_bin_countings,
                equal_var=False)
            output_fh.write("Margin: {}\n".format(self._margin))
            output_fh.write("Feature of interest median: {}\n".format(
                np.median(self._feature_bin_countings)))
            output_fh.write("Feature of interest mean: {}\n".format(
                np.mean(self._feature_bin_countings)))
            output_fh.write("Feature of interest standard deviation: "
                            "{}\n".format(np.std(self._feature_bin_countings)))
            output_fh.write("Background median: {}\n".format(
                np.median(self._random_bin_countings)))
            output_fh.write("Background mean: {}\n".format(
                np.mean(self._random_bin_countings)))
            output_fh.write("Background standard deviation: "
                            "{}\n".format(np.std(self._random_bin_countings)))
            output_fh.write("Welch's t-test statistic: {}\n".format(tstat))
            output_fh.write("Welch's t-test p-value: {}\n".format(pvalue))

    def write_countings_for_file(self):
        with open("{}_countings.txt".format(self._output_prefix),
                  "w") as output_fh:
            output_fh.write("Feature of interest ({}):\t{}\n".format(
                self._feature, ", ".join(
                    [str(value) for value in self._feature_bin_countings])))
            output_fh.write("Background:\t{}\n".format(", ".join(
                    [str(value) for value in self._random_bin_countings])))

    def plot_distribution(self):
        plt.style.use('ggplot')
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.hist(
            np.array(self._feature_bin_countings),
            alpha=0.5, bins=200,
            label="feature", normed=1)
        ax.hist(
            self._random_bin_countings, alpha=0.5, bins=200,
            label="background", normed=1)
        fig.tight_layout()
        ax.legend(loc='upper center')
        plt.savefig("{}_histograms.pdf".format(self._output_prefix))


def remove_position_information(name_with_pos_info):
    # Return just the chromosome part without the exact window
    # location
    return "-".join(name_with_pos_info.split("-")[:-1])

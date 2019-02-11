import gffutils
import numpy as np
import itertools
import matplotlib
matplotlib.use("Agg")
import hicsuntdracones.hicmatrix


class ColocalizationBase:

    def __init__(self, matrix_file, gff_file, feature, bin_size,
                 output_prefix, number_of_subsamplings, margin,
                 flanks_only, remove_zeros):
        self._matrix_file = matrix_file
        self._gff_file = gff_file
        self._feature = feature
        self._bin_size = int(bin_size)
        self._output_prefix = output_prefix
        self._number_of_subsamplings = int(number_of_subsamplings)
        self._margin = int(margin)
        self._flanks_only = flanks_only
        self._remove_zeros = remove_zeros
        self._features = None
        self._feature_overlapping_bins = []

    def read_gff_file(self):
        print("- Reading GFF file")
        self._features = gffutils.create_db(self._gff_file, ":memory:")

    def add_matrix_bins_as_features(self, matrix_file, features):
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
                end=pos_adjusted + self._bin_size,
                strand="+",
                attributes=f"ID={genome_bin}")])
        return interaction_matrix

    def extract_features_overlapping_bins(self):
        print("- Extracting bins that overlap given features")
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

    def _extract_feature_overlapping_bins_feature_and_flanks(self,
                                                             feature_start,
                                                             feature_end,
                                                             feature_seqid):
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

    def _extract_feature_overlapping_bins_flanks_only(self,
                                                      feature_start,
                                                      feature_end,
                                                      feature_seqid):
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
            interaction_countings = interaction_matrix.loc[
                [interaction_bin], other_bins]
            combined_interaction_countings.extend(
                list(itertools.chain(*interaction_countings.values.tolist())))
        combined_interaction_countings = np.array(
            combined_interaction_countings)
        return combined_interaction_countings

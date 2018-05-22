import gffutils
import numpy as np
import hicsuntdracones.hicmatrix


class Virtual4C(object):

    def __init__(self, matrix_file, gff_file, bin_size, track_name,
                 output_file):
        self._matrix_file = matrix_file
        self._gff_file = gff_file
        self._bin_size = bin_size
        self._track_name = track_name
        self._output_file = output_file

    def read_gff_file(self):
        print("- Reading GFF file")
        self._features = gffutils.create_db(self._gff_file, ":memory:")
        
    def read_matrix_file_and_add_as_features(self):
        """Read the HiC matrix file and add each bin as further feature in the
        annotation.
        """
        print("- Reading HiC matrix file")
        self.hic_matrix = hicsuntdracones.hicmatrix.HiCMatrix(
            self._matrix_file)
        self.hic_matrix.hic_matrix_df.set_index(
            self.hic_matrix.hic_matrix_df["Regions"], inplace=True)
        self.hic_matrix.normalize_by_columns_sum()
        for genome_bin in self.hic_matrix.bins():
            chrom = hicsuntdracones.hicmatrix.remove_position_information(
                genome_bin)
            # As the bin counting starts with 0 but the gff starts a 1
            # we have to add 1 here to the position.
            pos = hicsuntdracones.hicmatrix.bin_number(genome_bin)
            pos_adjusted = pos + 1
            self._features.update([gffutils.Feature(
                seqid=chrom,
                source='-',
                featuretype='HiC_bin',
                start=pos_adjusted,
                end=pos_adjusted+self._bin_size,
                strand="+",
                attributes="ID={}".format(genome_bin))])

    def extract_features_overlapping_bins(self):
        print("- Extracting bins that overlap given features")
        self._feature_overlapping_bins = []
        for feature_of_interest in self._features.all_features():
            if feature_of_interest.featuretype == "HiC_bin":
                continue
            for overlapping_bin in self._features.region(
                seqid=feature_of_interest.seqid,
                    start=feature_of_interest.start,
                    end=feature_of_interest.end,
                    featuretype="HiC_bin", completely_within=False):
                self._feature_overlapping_bins.append(
                    overlapping_bin.attributes["ID"][0])

    def generate_wiggle_file(self):
        """
        TODO
        - improve description
        - Use the bins of the feature of interest as selector of the
          series (column) and then the all the bin as selector of the
          row.

        """
        print("- Writing wiggle file")
        output_fh = open(self._output_file, "w")
        output_fh.write(
            'track type=wiggle_0 name="{}"\n'.format(self._track_name))
        current_chrom = ""
        for hic_bin_feature in self._features.all_features():
            if hic_bin_feature.featuretype != "HiC_bin":
                continue
            # In case a bin with a new chromosome shows up here, write
            # a new chrom definition line before writing the position
            # with the coverage value
            if hic_bin_feature.seqid != current_chrom:
                current_chrom = hic_bin_feature.seqid
                output_fh.write('variableStep chrom={} span={}\n'.format(
                    current_chrom, self._bin_size))
            hic_bin_name = hic_bin_feature.attributes["ID"][0]
            interaction_values = self.hic_matrix.matrix_values()[
                self._feature_overlapping_bins].loc[hic_bin_name].tolist()
            output_fh.write("{} {}\n".format(
                hic_bin_feature.start, np.mean(interaction_values)))
        output_fh.close()


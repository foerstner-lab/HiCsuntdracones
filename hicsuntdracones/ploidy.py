#!/usr/bin/env python
"""
FUNCTION: 

USAGE: 

Copyright (c) 2017, Konrad Foerstner <konrad@foerstner.org>

Permission to use, copy, modify, and/or distribute this software for
any purpose with or without fee is hereby granted, provided that the
above copyright notice and this permission notice appear in all
copies.

THE SOFTWARE IS PROVIDED 'AS IS' AND THE AUTHOR DISCLAIMS ALL
WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE
AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL
DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER
TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
PERFORMANCE OF THIS SOFTWARE.
         
"""
__description__ = ""
__author__ = "Konrad Foerstner <konrad@foerstner.org>"
__copyright__ = "2017 by Konrad Foerstner <konrad@foerstner.org>"
__license__ = "ISC license"
__email__ = "konrad@foerstner.org"
__version__ = ""

from collections import defaultdict
import pandas as pd


class Ploidy:

    def __init__(self,
                 hic_pro_matrix_values=None,
                 hic_pro_matrix_coordinates=None,
                 ploidy_file=None,
                 output_matrix_values=None,
                 output_matrix_coordinates=None,
                 bin_size=None):
        self._hic_pro_matrix_values = hic_pro_matrix_values
        self._hic_pro_matrix_coordinates = hic_pro_matrix_coordinates
        self._ploidy_file = ploidy_file
        self._output_matrix_values = output_matrix_values
        self._output_matrix_coordinates = output_matrix_coordinates
        self._bin_size = bin_size

    def main(self):
        factor_table = pd.read_table(self._ploidy_file)
        interaction_matrix = self._read_hic_pro_matrix(
            self._hic_pro_matrix_values,
            self._hic_pro_matrix_coordinates)
        interaction_matrix.set_index("Regions", inplace=True)
        ploidy_patches, ploidy_factors = self._contruct_bin_patches(
            interaction_matrix, factor_table, self._bin_size)
        for overlapping_bins, ploidy_factor in zip(ploidy_patches, ploidy_factors):
            interaction_matrix.loc[
                overlapping_bins, overlapping_bins] = interaction_matrix.loc[
                    overlapping_bins, overlapping_bins] * ploidy_factor
        self._write_matrix_in_hic_pro_format(
            interaction_matrix, self._bin_size,
            self._output_matrix_values, self._output_matrix_coordinates)

    def _contruct_bin_patches(self, interaction_matrix, factor_table, bin_size):
        """Look for overlaps of the given annotations with bins. A bin needs
        to overlap only partially with an annotated region.

              ---------------------------     Annotation
            ====           ====         ====  Bin accepted as overlapping

        """

        # Build matrix of bin position to make querying easy
        bin_pos_matrix = pd.DataFrame()
        bin_pos_matrix.set_index = interaction_matrix.index
        bin_pos_matrix["HiCMatrix"] = interaction_matrix["HiCMatrix"]
        bin_pos_matrix["chrom"] = bin_pos_matrix["HiCMatrix"].apply(
            lambda bin_name: "-".join(bin_name.split("-")[:-1]))
        bin_pos_matrix["start"] = bin_pos_matrix["HiCMatrix"].apply(
            lambda bin_name: int(bin_name.split("-")[-1]) + 1)
        bin_pos_matrix["end"] = bin_pos_matrix["start"] + bin_size - 1
        bin_pos_matrix["start"] = bin_pos_matrix["start"].astype(int)
        bin_pos_matrix["end"] = bin_pos_matrix["end"].astype(int)

        ploidy_patches = []  # Will become a list of lists
        ploidy_factors = factor_table["ploidy_factor"].tolist()
        for index, patch in factor_table.iterrows():
            overlapping_bins = bin_pos_matrix[
                (bin_pos_matrix["chrom"] == patch["chrom"])
                &
                (
                    (
                        (bin_pos_matrix["start"] >= patch["start"])
                        &
                        (bin_pos_matrix["start"] <= patch["stop"])
                    )
                    |
                    (
                        (bin_pos_matrix["end"] >= patch["start"])
                        &
                        (bin_pos_matrix["end"] <= patch["stop"])
                    )
                )
            ]
            ploidy_patches.append(overlapping_bins["HiCMatrix"].tolist())
        assert len(ploidy_patches) == len(ploidy_factors)
        return(ploidy_patches, ploidy_factors)

    def _write_matrix_in_hic_pro_format(self,
                                        interaction_matrix,
                                        bin_size,
                                        output_matrix_coordinates,
                                        output_matrix_values):

        coordinates = pd.DataFrame()
        coordinates["Regions"] = interaction_matrix.index
        coordinates["Chrom_name"] = coordinates["Regions"].apply(
            lambda bin_name: "-".join(bin_name.split("-")[:-1]))
        coordinates["Start_pos"] = coordinates["Regions"].apply(
            lambda bin_name: int(bin_name.split("-")[-1]))
        coordinates["End_pos"] = coordinates["Start_pos"].apply(
            lambda start_pos: start_pos + bin_size)
        coordinates["Id"] = [start_pos + 1 for start_pos in range(
            len(coordinates["Regions"]))]
        del coordinates["Regions"]
        coordinates.to_csv(
            output_matrix_coordinates, index=False, sep="\t", header=False)

        # Generate non redundant comparison list by going throug the matrix
        # like this:
        # |
        # ||
        # |||
        # ||||
        # |||||
        # ||||||
        del interaction_matrix["HiCMatrix"]
        start_row_number = 0
        with open(output_matrix_values, "w") as output_fh:
            for column_number, bin_name in enumerate(interaction_matrix.columns):
                for row_number in range(
                        start_row_number, len(interaction_matrix.columns)):
                    interaction_value = int(round(interaction_matrix.iloc[
                        row_number, column_number]))
                    if interaction_value == 0:
                        continue
                    output_fh.write("{}\t{}\t{}\n".format(
                        column_number + 1, row_number + 1,
                        interaction_value))
                start_row_number += 1

    def _read_hic_pro_matrix(self, matrix_values_file, matrix_coordinates_file):
        pair_value_table = pd.read_table(
            matrix_values_file, names=["bin_a", "bin_b", "counting"])
        binning_information = pd.read_table(
            matrix_coordinates_file, names=["replicon", "start", "end", "bin_id"])

        bin_id_to_name = dict([
            (bin_id, "-".join([replicon, str(start)]))
            for bin_id, replicon, start in zip(
                    binning_information.bin_id,
                    binning_information.replicon,
                    binning_information.start)])

        bin_pair_to_value = defaultdict(dict)
        for bin_a, bin_b, value in zip(
                pair_value_table["bin_a"], pair_value_table["bin_b"],
                pair_value_table["counting"]):
            bin_pair_to_value[bin_a][bin_b] = value
            bin_pair_to_value[bin_b][bin_a] = value
        result_matrix = pd.DataFrame()
        result_matrix["HiCMatrix"] = binning_information["bin_id"].apply(
            lambda bin_id: bin_id_to_name[bin_id])
        result_matrix["Regions"] = result_matrix["HiCMatrix"]

        # Add column for each bin
        for bin_a in binning_information["bin_id"]:
            result_matrix[bin_id_to_name[
                bin_a]] = binning_information["bin_id"].apply(
                lambda bin_b: bin_pair_to_value[bin_a].get(bin_b, 0))
        return result_matrix


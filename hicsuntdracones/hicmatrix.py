import sys
import numpy as np
import pandas as pd


class HiCMatrix():

    def __init__(self, hic_matrix_file: str):
        self._hic_matrix_file = hic_matrix_file
        self._read_matrix()

    def _read_matrix(self):
        self._check_file()
        self.hic_matrix_df = pd.read_table(self._hic_matrix_file)
        self._check_hic_matrix_df()
        self.number_of_bins = self.hic_matrix_df.shape[0]
        self.chromosomes = self._chromosomes()
        
    def _check_file(self):
        # TODO
        # check if file contains "HiCMatrix" and "Regions"
        pass

    def _check_hic_matrix_df(self):
        no_of_rows, no_of_columns = self.hic_matrix_df.shape
        if not no_of_columns - 2 == no_of_rows:
            sys.stderr.write("Unexpected ratio of columns and rows."
                             "Is this really a HiC matrix in Homer format?")
            sys.exit(1)
    
    def save(self, output_hic_matrix_file):
        self.hic_matrix_df.to_csv(output_hic_matrix_file, sep="\t",
                                  index=False)

    def normalize_by_columns_sum(self):
        """
        Assumption: input is an iced matrix.
        """
        column_sum_median = self._calc_column_sum_median()
        for column in self.hic_matrix_df.columns:
            if column in ["HiCMatrix", "Regions"]:
                continue
            self.hic_matrix_df[column] = self.hic_matrix_df[
                column] / column_sum_median

    def _calc_column_sum_median(self):
        """Needed as the sums are not completely identical and some row sums
        are 0. Assumption: input is an iced matrix.
        """
        return np.median([
            self.hic_matrix_df[col].sum()
            for col in self.hic_matrix_df.columns[2:]])

    def _chromosomes(self):
        return sorted(self.hic_matrix_df[
            "Regions"].apply(remove_position_information).unique())


def remove_position_information(name_with_pos_info):
    # Return just the chromosome part without the exact window
    # location
    return "-".join(name_with_pos_info.split("-")[:-1])


import sys
import pandas as pd


class HiCMatrix():

    def __init__(self, hic_matrix_file):
        self._hic_matrix_file = hic_matrix_file
        self._read_matrix()

    def _read_matrix(self):
        self._check_file()
        self.hic_matrix_df = pd.read_table(self._hic_matrix_file)
        self._check_hic_matrix_df()
        self.hic_matrix_df["Regions"] = self.hic_matrix_df[
            "Regions"].apply(remove_position_information)
        self.number_of_bins = self.hic_matrix_df.shape[0]
        
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
    
    def to_csv():
        pass

    def _chromosomes():
        pass


def remove_position_information(name_with_pos_info):
    # Return just the chromosome part without the exact window
    # location
    return "-".join(name_with_pos_info.split("-")[:-1])
    

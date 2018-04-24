import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import ndimage


class HiCMatrix():

    def __init__(self, hic_matrix_file=None, hic_matrix_df=None):
        if hic_matrix_file is not None:
            self._hic_matrix_file = hic_matrix_file
            self._read_matrix()
        elif hic_matrix_df is not None:
            self.hic_matrix_df = hic_matrix_df

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
        if not self.hic_matrix_df.columns[0] == "HiCMatrix":
            sys.stderr.write("Missing column 'HiCMatrix'."
                             "Is this really a HiC matrix in Homer format?")
            sys.exit(1)
        if not self.hic_matrix_df.columns[1] == "Regions":
            sys.stderr.write("Missing column 'Regions'."
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

    def select(self, keep_pattern=None, remove_pattern=None, inplace=False):
        """Select a submatrix based on given filters of the bin names.

        We like minimalism => Removing is stronger than keeping.
        """
        if inplace:
            self.hic_matrix_df = self._select(
                keep_pattern=keep_pattern, remove_pattern=remove_pattern)
        else:
            return HiCMatrix(hic_matrix_df=self._select(
                keep_pattern=keep_pattern, remove_pattern=remove_pattern))

    def _select(self, keep_pattern=None, remove_pattern=None):
        filtered_bins = self.bins()
        if keep_pattern is not None:
            filtered_bins = filtered_bins[
                filtered_bins.str.contains(keep_pattern)]
        if remove_pattern is not None:
            filtered_bins = filtered_bins[
                ~ filtered_bins.str.contains(remove_pattern)]
        # Filter columns
        submatrix = self.hic_matrix_df[
            ["HiCMatrix", "Regions"] + filtered_bins.tolist()]
        # Filter rows
        submatrix = submatrix[
            submatrix["HiCMatrix"].isin(filtered_bins.tolist())]
        #  As the index still still based on the odering before
        #  removing row the index has to be written.
        submatrix.index = range(len(submatrix))
        return submatrix

    def bins(self):
        return self.hic_matrix_df["Regions"].rename("bins", inplace=True)

    def div_by(self, denominator_matrix, pseudocount=0.001, inplace=False):
        if inplace:
            self.hic_matrix_df = self._div_by(
                denominator_matrix, pseudocount=pseudocount)
        else:
            return HiCMatrix(hic_matrix_df=self._div_by(
                denominator_matrix, pseudocount=pseudocount))

    def _div_by(self, denominator_matrix, pseudocount) -> pd.DataFrame:
        """

        Add pseudocount first and then normalize to make sure that
        column sums are the same.
        """
        numerator_matrix_values = self.matrix_values() + pseudocount
        denominator_matrix_values = (
            denominator_matrix.matrix_values() + pseudocount)
        diff_matrix = numerator_matrix_values / denominator_matrix_values
        return pd.concat([self.hic_matrix_df[[
            "HiCMatrix", "Regions"]], diff_matrix],
                         axis=1, join_axes=[self.hic_matrix_df.index])

    def matrix_values(self):
        """Return the matrix without the bin name columns.
        """
        return self.hic_matrix_df.ix[0:, 2:]

    def heatmap(self, vmin=None, vmax=None, by_chrom=False, rotate=False,
                output_pdf=None, output_png_prefix=None, png_dpi=600):
        self._vmax = vmax
        self._vmin = vmin
        self._output_pdf = output_pdf
        self._output_png_prefix = output_png_prefix
        self._pp = None
        self._png_dpi = png_dpi
        self._rotate = rotate
        if not by_chrom:
            self._plot_heatmap_globally()
        else:
            self._plot_heatmap_split_by_chrom()
        if self._output_pdf is not None:
            self._pp.close()

    def _plot_heatmap_globally(self):
        self._plot_heatmap(self)

    def _plot_heatmap_split_by_chrom(self):
        chroms = sorted(set([
            "-".join(genome_bin.split("-")[:-1])
            for genome_bin in self.matrix_values().columns]))
        for chrom in chroms:
            chrom_hic_matrix = self.select(keep_pattern=chrom)
            # Make sure that at lest 2 bin are in the submatrix
            if chrom_hic_matrix.hic_matrix_df.shape[0] < 2:
                print("Skipping {} as the number of bins is too low.".format(
                    chrom))
                continue
            self._plot_heatmap(chrom_hic_matrix, title=chrom)
    
    def _plot_heatmap(self, hic_matrix, title=""):
        matrix_values = hic_matrix.matrix_values()
        if self._rotate:
            matrix_values = ndimage.rotate(matrix_values, 45.0)
            matrix_values = matrix_values[
                :int(matrix_values.shape[0]/2), :]
        heatmap = sns.heatmap(
            matrix_values, square=True,
            vmax=self._vmax, vmin=self._vmin,
            xticklabels=False, yticklabels=False,
            cmap=sns.cubehelix_palette(n_colors=500))
        plt.title(title)
        plt.tight_layout()
        self._write_heatmap_to_file(heatmap, title=title)
        plt.close()

    def _write_heatmap_to_file(self, heatmap, title=""):
        if self._output_pdf is not None:
            if self._pp is None:
                self._pp = PdfPages(self._output_pdf)
            self._pp.savefig(heatmap.figure)
        elif self._output_png_prefix is not None:
            if title == "":
                output_file_name = self._output_png_prefix + ".png"
            else:
                output_file_name = "{}_{}.png".format(
                    self._output_png_prefix, title)
            heatmap.figure.savefig(output_file_name, dpi=self._png_dpi)
        else:
            raise MissingOutputFile


def remove_position_information(name_with_pos_info: str):
    # Return just the chromosome part without the exact window
    # location
    return "-".join(name_with_pos_info.split("-")[:-1])


def read_hic_matrix(input_file: str):
    return HiCMatrix(hic_matrix_file=input_file)


def MissingOutputFile(BaseException):
    pass

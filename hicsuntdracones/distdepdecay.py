from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import hicsuntdracones.hicmatrix


class DistDepDecayOutputGenerator():

    def __init__(self, matrix_file, bin_size, output_prefix):
        hic_matrix = hicsuntdracones.hicmatrix.HiCMatrix(matrix_file)
        self._output_prefix = output_prefix
        self._chroms_dists_and_countings = (
            hic_matrix.calc_distance_dependent_decay(bin_size))
        # Just a different arrangement of the data:
        self._chroms_countings_by_dist = defaultdict(
            lambda: defaultdict(list))
        for chrom in self._chroms_dists_and_countings.keys():
            for dist, counting in zip(
                    self._chroms_dists_and_countings[chrom]["dists"],
                    self._chroms_dists_and_countings[chrom]["countings"]):
                self._chroms_countings_by_dist[chrom][dist].append(
                    counting)

    def write_table_file(self):
        with open("{}_distance_dependent_decay.csv".format(
                self._output_prefix), "w") as output_fh:
            output_fh.write("Chrom\tDist\tMean counting value\t"
                            "Counting standard deviation\t"
                            "Countings\tNumber of countings\n")
            for chrom in self._chroms_dists_and_countings.keys():
                for dist in sorted(self._chroms_countings_by_dist[
                        chrom].keys()):
                    countings = self._chroms_countings_by_dist[chrom][dist]
                    counting_mean = np.mean(countings)
                    counting_stand_dev = np.std(countings)
                    output_fh.write(
                        "{}\t{}\t{}\t{}\t{}\t{}\n".format(
                            chrom, dist, counting_mean,
                            counting_stand_dev,
                            ", ".join([str(count) for count in countings]),
                            len(countings)))

    def plot_bin_averages(self):
        _pp_av = PdfPages(
            "{}_distance_dependent_decay_averages.pdf".format(
                self._output_prefix))
        plt.style.use('ggplot')
        for chrom in self._chroms_countings_by_dist.keys():
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            dists = sorted(self._chroms_countings_by_dist[chrom].keys())
            counting_means = [np.mean(self._chroms_countings_by_dist[
                chrom][dist]) for dist in dists]
            ax.plot(dists, counting_means, "-", alpha=0.5, color="black")
            # ax.set_ylim([0, 20])
            ax.set_title(chrom)
            _pp_av.savefig(fig)
            plt.close(fig)
        _pp_av.close()

    def plot_bin_averages_with_error_bars(self):
        _pp_av = PdfPages(
            "{}_distance_dependent_decay_averages_with_error_bars.pdf".format(
                self._output_prefix))
        plt.style.use('ggplot')
        for chrom in self._chroms_countings_by_dist.keys():
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            dists = sorted(self._chroms_countings_by_dist[
                chrom].keys())
            counting_means = [np.mean(self._chroms_countings_by_dist[
                chrom][dist]) for dist in dists]
            counting_stand_dev = [np.std(self._chroms_countings_by_dist[
                chrom][dist]) for dist in dists]
            ax.errorbar(dists, counting_means, yerr=counting_stand_dev,
                        color="black")
            # ax.set_xscale("log")
            # Might raise error with log values
            try:
                ax.set_yscale("log")
                # ax.set_ylim([0, 20])
                ax.set_title(chrom)
                _pp_av.savefig(fig)
                plt.close(fig)
            except ValueError:
                pass
            
        _pp_av.close()
                    
    def plot_all_values(self):
        self._pp = PdfPages("{}_distance_dependent_decay_"
                            "all_values.pdf".format(self._output_prefix))
        plt.style.use('ggplot')
        for chrom in self._chroms_dists_and_countings.keys():
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax.plot(
                self._chroms_dists_and_countings[chrom]["dists"],
                self._chroms_dists_and_countings[chrom]["countings"],
                ".", alpha=0.5)
            ax.set_title(chrom)
            self._pp.savefig(fig)
            plt.close(fig)
        self._pp.close()
        


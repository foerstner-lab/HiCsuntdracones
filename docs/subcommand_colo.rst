====
colo
====

The colocalization tester allows to select regions of interest and to
test them for increased interchromosomal interaction frequencies in
two ways.  Either, genomic loci of interest can be selected by
keywords/features from the corresponding genome gff file (._feature),
which are then translated into the corresponding bins of the
interaction matrix and collected. Alternatively, genomic regions of
any length and choice can be subjected to the colocalization tester in
gff format (.read_gff_file), and overlapping bins will be compiled
from the matrix (.read_matrix_file_and_add_as_features) in the same
way. The regions of interest will be tested against a randomized set
of bins (.extract_random_bins), that mirrors the query sample in
distribution and length of spanned bins. The distributions of query
and randomized sample are plotted (.plot_distribution) and tested for
significant differences (.perform_t_test) using Welch’s t-test
statistics. This sampling can be redone as many times as desired
(._number_of_subsamplings) and bins with zero values can be in- or
excluded from the analysis (._remove zeros). Additionally, the tool
allows to add various numbers of flanking bins to the regions of
interest (._margin) and also to restrict analysis to these
(._flanks_only). As output, the tool gives the individual interaction
counts from the query and the randomized sample, a histogram of those
as a pdf, as well as mean and median values with according standard
deviation and p-values.

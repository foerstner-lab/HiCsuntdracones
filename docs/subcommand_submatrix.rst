=========
submatrix
=========

HiCMatrix first tests the shape of the converted matrix. To be able to
compare different iced matrices (with different column sums),
normalize_by_columns_sum allows to divide each column by the column
sum median, thereby turning the interaction values into
frequencies. In some cases, it may be useful to restrict the analysis
to a subset of the matrix (select). Thereby, distinct chromosomes can
be in- (keep_pattern) or excluded (remove_pattern) from the overall
matrix, resulting in a submatrix that can be saved separately or used
instead of the original matrix for the following analyses (inplace).
Hi-C Matrix also gives the possibility to generate differential
matrices by dividing the matrix (or submatrix) by a denominator
matrix. As some bins contain zero values, a pseudocount of choice can
be added to each bin of the nominator and denominator matrix (default:
0.001).

::

    $ hicsd submatrix -h

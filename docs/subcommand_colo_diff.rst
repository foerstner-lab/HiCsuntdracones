=========
colo_diff
=========

Differential colocalization can be used to compare distinct
interchromosomal interactions in two different matrices. It assesses
the interaction frequencies of genomic loci of interest over a
randomized background sample (as described in Colocalization) in two
matrices of choice and tests for significant differences in
signal-over-background among them. In each round of subsampling, the
query sample is tested against a randomized background sample in both
matrices. The ratio between average signal and average background is
compiled and the subsampling is repeated (as often as desired), thus
resulting in two respective populations of signal-over-background
ratios (one for each matrix). The difference between these two
populations is tested using Welch’s t-test statistics.

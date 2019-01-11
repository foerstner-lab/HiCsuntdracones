==========
virtual_4C
==========

This tool enables to use any region of choice from the full-genome HiC
matrix as bait for a virtual 4C. The regions of choice can be
submitted in gff format (see example: Chromosome – start –
end). Subsequently, these genomic regions are translated into the
corresponding bins of the HiC matrix and the affected columns are
extracted. If several columns were selected, an average interaction
value for every genomic bin is calculated. The interaction values for
each bin are translated into a wiggle file, where step and window size
correspond to the resolution of the input matrix.

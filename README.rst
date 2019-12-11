=======================================================
HiC sunt dracones - Your little helper for HiC analyses
=======================================================


.. image:: https://img.shields.io/pypi/v/hicsuntdracones.svg
        :target: https://pypi.python.org/pypi/hicsuntdracones

.. image:: https://img.shields.io/travis/foerstner-lab/hicsuntdracones.svg
        :target: https://travis-ci.org/foerstner-lab/hicsuntdracones

.. image:: https://readthedocs.org/projects/hicsuntdracones/badge/?version=latest
        :target: https://hicsuntdracones.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://pyup.io/repos/github/foerstner-lab/hicsuntdracones/shield.svg
     :target: https://pyup.io/repos/github/foerstner-lab/hicsuntdracones/
     :alt: Updates

.. image:: https://zenodo.org/badge/95952483.svg
     :target: https://zenodo.org/badge/latestdoi/95952483
     :alt: Zenodo
	   
.. image:: https://upload.wikimedia.org/wikipedia/commons/2/27/1601_De_Bry_and_de_Veer_Map_of_Nova_Zembla_and_the_Northeast_Passage_-_Geographicus_-_NovaZembla-debry-1601.jpg
   :height: 75px

Welcome to HiC sunt dracones (or HiCSD in short).
	    
-------------------------------------
Why should you use HiC sunt dracones?
-------------------------------------

The functions of HiCsd are built on both Hi-C Pro matrices and Homer matrices. HiCsd uses HiCUP for truncation and filtering of reads and implicit matrix normalization (ICEing) as implemented Hi-C Pro. It is further possible, to explicitly normalize read counts of certain regions in the raw matrix for relative differences in ploidy. HiCsd then allows to analyse the full HiC matrix, differential matrices, a desired subset of a matrix (submatrix) as well as individual chromosomes. These matrices can be visualized either as full heatmaps or triangles with color and scale of the heatmaps adjustable. Using HiCsd, testing colocalization of loci of interest is easy as HiCsd can pick these loci from a genome annotation by searching for given keywords or by feeding a list of loci of interest into the program. Also, HiCsd is able to compare the interaction frequencies of distinct loci between different matrices. 
The distance dependent decay of interaction frequencies can be calculated for both, individual chromosomes or as a mean/median function including all chromosomes of choice. Also, HiCsd is able to extract the interaction frequencies between any locus of choice with all other genomic loci, thus generating 4C-like profiles in wiggle format, that can easily be explored and compared using a common genome browser.

	    
-------------------
Subcommand overview
-------------------

::
    
    usage: hicsd [-h] {version} ...
    
    positional arguments:
      {version}   commands
        version   Show version
    
    optional arguments:
      -h, --help  show this help message and exit


    usage: hicsd [-h]
                 {version,hicpro2homer,number_of_bins,chromosomes,norm_by_col_sum,submatrix,diff_matrix,heatmap,virtual_4C,dist_dep_decay,colo,colo_diff}
                 ...

    positional arguments:
      {version,hicpro2homer,number_of_bins,chromosomes,norm_by_col_sum,submatrix,diff_matrix,heatmap,virtual_4C,dist_dep_decay,colo,colo_diff}
                            commands
        version             Show version
        hicpro2homer        Convert a matrix in HiC-Pro format to Homer format
        number_of_bins      Return the number of bins of the matrix
        chromosomes         Return the chromosomes used in the matrix
        norm_by_col_sum     Normalize by column sum to make matrices comparable
        submatrix           Extract submatrices
        diff_matrix         Generate differential matrix by dividing the values of
                            one matrix by the values of the other.
        heatmap             Plot interaction matix heatmap
        virtual_4C          Perform virtual 4C analysis
        dist_dep_decay      Distant dependent decay
        colo                Perform colocalisation analysis
        colo_diff           Perform differential colocalisation analysis

    optional arguments:
      -h, --help            show this help message and exit
      

-----------
Development
-----------

This repo contains 2 branches: "master" for development and "production" for stable releases


------
Trivia
------

- `Origin of the name <https://en.wikipedia.org/wiki/Here_be_dragons>`__

=======================================================
HiC sunt dracones - Your little helper for HiC analyses
=======================================================


.. image:: https://img.shields.io/pypi/v/hicsuntdracones.svg
        :target: https://pypi.python.org/pypi/hicsuntdracones

.. image:: https://img.shields.io/travis/konrad/hicsuntdracones.svg
        :target: https://travis-ci.org/konrad/hicsuntdracones

.. image:: https://readthedocs.org/projects/hicsuntdracones/badge/?version=latest
        :target: https://hicsuntdracones.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://pyup.io/repos/github/konrad/hicsuntdracones/shield.svg
     :target: https://pyup.io/repos/github/konrad/hicsuntdracones/
     :alt: Updates

.. image:: https://upload.wikimedia.org/wikipedia/commons/2/27/1601_De_Bry_and_de_Veer_Map_of_Nova_Zembla_and_the_Northeast_Passage_-_Geographicus_-_NovaZembla-debry-1601.jpg
   :height: 75px

-------------------------------------
Why should you use HiC sunt dracones?
-------------------------------------


	    
-----
Usage
-----

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
      
------
Trivia
------

- `Origin of the name <https://en.wikipedia.org/wiki/Here_be_dragons>`__

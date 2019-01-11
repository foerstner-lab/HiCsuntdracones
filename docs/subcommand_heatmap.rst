=======
heatmap
=======

The subcommand heatmap prints matrices, submatrices or differential
matrices of choice, either as whole genome heatmaps
(_plot_heatmap_globally) or chromosome-wise
(_plot_heatmap_split_by_chrom; two bins required as minimum
lenght). The heatmap can either be plotted in a square format or as
triangles (--rotate). Colors are based on the colour palettes of
seaborn and can be adjusted as desired (default: sns.cubehelix_palette
for non-differential heatmaps and Rdblue for differential heatmaps)
and the minimum and maximum values of the color scale can be adjusted
(._vim, ._vmax). The output file (_write_heatmap_to_file) can either
be in png or pdf format. In png format, all generated heatmaps are
printed in separate files, whereas in pdf format, all files are
printed to a single file. Also, the png file resolution can be set as
desired (._png_dpi).


::

     $ heatmap -h
     usage: hicsd heatmap [-h] --matrix MATRIX_FILE [--vmin VMIN] [--vmax VMAX]
                          [--rotate] [--by_chrom] [--output_pdf OUTPUT_PDF]
                          [--output_png_prefix OUTPUT_PNG_PREFIX]
                          [--png_dpi PNG_DPI] [--differential]
     
     optional arguments:
       -h, --help            show this help message and exit
       --matrix MATRIX_FILE, -m MATRIX_FILE
       --vmin VMIN, -i VMIN
       --vmax VMAX, -a VMAX
       --rotate, -r          Rotate 45 degrees. Pyramids FTW!
       --by_chrom, -c
       --output_pdf OUTPUT_PDF, -p OUTPUT_PDF
                             Name of the output file. All figures will be written
                             into this single file
       --output_png_prefix OUTPUT_PNG_PREFIX, -n OUTPUT_PNG_PREFIX
                             Prefix for output png files. Each figure file be
                             stored in a separate file.
       --png_dpi PNG_DPI, -e PNG_DPI
       --differential, -d

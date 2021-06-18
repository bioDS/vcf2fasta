# Run VCF files with Beast
Transform VCF files into fasta and run them using the NucleotideDiploid16 model in BEAST

## Requirement:
**Python 3** with `pysam` and `snakemake` packages. These can be installed with `python3 -m pip install pysam snakemake`.
 
**R** with beter, phyloRNA, argparser. Use `devtools` to install the `beter` and `phyloRNA` packages `devtools::install_github(c("biods/beter", "biods/phyloRNA"))`. The `argparser` package can be installed from cran `install.packages("argparser")`.
 
**BEAST2** and the `phylonco` package. To install `phylonco` using the `BEAST2` package manager, type: `packagemanager -add phylonco` or use the `beauti` GUI.

## Run the analysis
To run the analysis, navigate to the directory with the `snakefile` and type `snakemake`. In some versions of the `snakemake`, you might have to specify the number of CPU cores used during the analysis, such as: `snakemake -j 4` for 4 CPU cores.

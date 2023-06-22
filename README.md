# Run VCF files with Beast
Transform VCF files into nexus or fasta format for use with the BEAST phylonco package.

## Requirement:
**File conversion:**

* Python 3 with `pysam` and `snakemake` packages. These can be installed with `python3 -m pip install pysam snakemake`.
 
**Creating BEAST XMLs:**

* Beauti GUI can be used to generate XMLs

* Alternatively, for BEAST 2.6 XMLs can be generated using R with [beter](https://github.com/bioDS/beter/tree/master), [phyloRNA](https://github.com/bioDS/phyloRNA), argparser. Use `devtools` to install the `beter` and `phyloRNA` packages `devtools::install_github(c("biods/beter", "biods/phyloRNA"))`. The `argparser` package can be installed from cran `install.packages("argparser")`. 
 
**BEAST requirements:** 

* BEAST 2 and the `phylonco` package. To install `phylonco` using the `BEAST2` package manager, type: `packagemanager -add phylonco` or use the `beauti` GUI.

## Example VCF files
Example VCF files can be downloaded from https://github.com/amkozlov/cellphy-paper/blob/main/data/

The examples below use the VCF file `E15.CellPhy.vcf` from [E15.zip](https://github.com/amkozlov/cellphy-paper/blob/main/data/E15.zip)

Citation for dataset: https://doi.org/10.1186/s13059-021-02583-w

## Converting VCF to nexus

**Full diploid genotype**

Replace `<input vcf file>` with the path of your VCF, and `<output nexus file>` with the path of your Nexus file. 

```
python3 src/vcf2nexus.py <input vcf file> <output nexus file> --encoding nd16
```

Example usage:
```
python3 src/vcf2nexus.py E15.CellPhy.vcf E15.CellPhy.nex --encoding nd16
```

**Binary genotype**
```
python3 src/vcf2nexus.py <input vcf file> <output nexus file> --encoding binary
```

Example usage:
```
python3 src/vcf2nexus.py E15.CellPhy.vcf E15.CellPhy.binary.nex --encoding binary
```

## Converting VCF to fasta

**Full diploid genotype**

Replace `<input vcf file>` with the path of your VCF, and `<output fasta file>` with the path of your Fasta file. 

```
python3 src/vcf2fasta.py <input vcf file> <output fasta file> --encoding nd16
```

Example usage:
```
python3 src/vcf2fasta.py E15.CellPhy.vcf E15.CellPhy.fasta --encoding nd16
```

**Binary genotype**
```
python3 src/vcf2fasta.py <input vcf file> <output fasta file> --encoding binary
```

Example usage:
```
python3 src/vcf2fasta.py E15.CellPhy.vcf E15.CellPhy.binary.fasta --encoding binary
```

## Running analyses 
To run analyses using snakemake, navigate to the directory with the `snakefile` and type `snakemake`. In some versions of the `snakemake`, you might have to specify the number of CPU cores used during the analysis, such as: `snakemake -j 4` for 4 CPU cores.

To prevent `snakemake` deleting output of failed runs (which can be quite costly), run `snakemake` with the `--keep-incomplete` switch, e.g.:
```
snakemake -j 4 --keep-incomplete
```

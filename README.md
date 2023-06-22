# VCF to Nexus or FASTA 
Convert VCF files into nexus or fasta formats for use with the BEAST phylonco package https://github.com/bioDS/beast-phylonco

## Requirements:
Python 3 with the `pysam` package. This can be installed with `python3 -m pip install pysam`.
 
## Example VCF files
Example VCF files can be downloaded from https://github.com/amkozlov/cellphy-paper/blob/main/data/

The examples below use the VCF file `E15.CellPhy.vcf` from [E15.zip](https://github.com/amkozlov/cellphy-paper/blob/main/data/E15.zip)

Citation for dataset: https://doi.org/10.1186/s13059-021-02583-w

## Converting VCF to nexus

**Full diploid genotype**

Replace `<vcf file>` with the path of your VCF, and `<nexus file>` with the path of your Nexus file. 

```
python3 src/vcf2nexus.py <vcf file> <nexus file> --encoding nd16
```

Example usage:
```
python3 src/vcf2nexus.py E15.CellPhy.vcf E15.CellPhy.nex --encoding nd16
```

**Binary genotype**
```
python3 src/vcf2nexus.py <vcf file> <nexus file> --encoding binary
```

Example usage:
```
python3 src/vcf2nexus.py E15.CellPhy.vcf E15.CellPhy.binary.nex --encoding binary
```

## Converting VCF to fasta

**Full diploid genotype**

Replace `<vcf file>` with the path of your VCF, and `<fasta file>` with the path of your Fasta file. 

```
python3 src/vcf2fasta.py <vcf file> <fasta file> --encoding nd16
```

Example usage:
```
python3 src/vcf2fasta.py E15.CellPhy.vcf E15.CellPhy.fasta --encoding nd16
```

**Binary genotype**
```
python3 src/vcf2fasta.py <vcf file> <fasta file> --encoding binary
```

Example usage:
```
python3 src/vcf2fasta.py E15.CellPhy.vcf E15.CellPhy.binary.fasta --encoding binary
```

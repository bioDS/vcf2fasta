#'/usr/bin/env python3
"""Parse a multi-sample VCF file and produce a FASTA file."""

import textwrap
import argparse
import pysam


nd16phased = {
    "AA" : "0",
    "AC" : "1",
    "AG" : "2",
    "AT" : "3",
    "CA" : "4",
    "CC" : "5",
    "CG" : "6",
    "CT" : "7",
    "GA" : "8",
    "GC" : "9",
    "GG" : "a",
    "GT" : "b",
    "TA" : "c",
    "TC" : "d",
    "TG" : "e",
    "TT" : "f",
    # Phased with one unknown
    ".A" : "g",
    ".C" : "h",
    ".G" : "i",
    ".T" : "j",
    "A." : "n",
    "C." : "o",
    "G." : "p",
    "T." : "q",
    ".." : "?"
    }


nd16unphased = {
    "AA" : "0",
    "CC" : "5",
    "GG" : "a",
    "TT" : "f",
    "AC" : "M",
    "CA" : "M",
    "AG" : "R",
    "GA" : "R",
    "AT" : "W",
    "TA" : "W",
    "CG" : "S",
    "GC" : "S",
    "CT" : "Y",
    "TC" : "Y",
    "GT" : "K",
    "TG" : "K",
    ".." : "?"
    }


def parse_args():
    parser = argparse.ArgumentParser("Transform a multi-sample VCF file into a Fasta format")
    parser.add_argument("vcf", type=str, help="A multi-sample VCF file")
    parser.add_argument("fasta", type=str, help="A path to an output fasta file")
    args = parser.parse_args()
    return args


def vcf2fasta(vcf, fasta):
    names, sequences = parse_vcf(vcf)
    write_fasta(names, sequences, fasta)


def write_fasta(names, sequences, file):
    items = ["\n".join([">" + name, textwrap.fill(seq, width=80)]) for name, seq in zip(names, sequences)]
    with open(file, "w") as f:
        f.write("\n\n".join(items))
        f.write("\n")

# TODO: Wrapper class for vcf would make this easier
def parse_vcf(file):
    with pysam.VariantFile(file) as vcf:
        names = list(vcf.header.samples)
        sequences = get_sequences(vcf.fetch())
    return names, sequences


def get_sequences(variants):
    sequences = map(variant2bases, variants)
    sequences = transpose_list(sequences)
    sequences = map("".join, sequences)
    return list(sequences)


def transpose_list(lst):
    return list(map(list, zip(*lst)))

# TODO: Wrapper class for variant would make this easier
def variant2bases(variant):
    return [sample2base(sample) for sample in variant.samples.itervalues()]

# TODO: Wrapper class for samples would make this easier
def sample2base(sample):
    return translate_alleles(sample.alleles, sample.phased)


def translate_alleles(alleles, phased=True):
    genome = ["." if allele is None else allele for allele in alleles]
    return translate_genome("".join(genome), phased=phased)


def translate_genome(genome, phased=True):
    if len(genome) != 2:
        raise ValueError(f"Genome must be diploid. The ploidy is: {len(genome)}")
    if phased:
        return nd16phased[genome]
    else:
        return nd16unphased[genome]


if __name__ == "__main__":
    args = parse_args()
    vcf2fasta(args.vcf, args.fasta)

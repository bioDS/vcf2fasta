#'/usr/bin/env python3
"""Parse a multi-sample VCF file and produce a Nexus file."""

import textwrap
import argparse
import pysam


datatype = {
    "binary": "binary",
    "nd16": "nucleotideDiploid16"
    }


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


binary = {
    "00" : "0",
    "11" : "1",
    "01" : "1",
    "10" : "1",
    ".1" : "1",
    "1." : "1",
    "0." : "?",
    ".0" : "?",
    ".." : "?"
    }


def parse_args():
    parser = argparse.ArgumentParser("Transform a multi-sample VCF file into a Nexus format")
    parser.add_argument("vcf", type=str, help="A multi-sample VCF file")
    parser.add_argument("nexus", type=str, help="A path to an output nexus file")
    parser.add_argument("--encoding", type=str, choices=["nd16", "binary"], default="nd16",
        help="""
    A type of encoding used during to translate a diploid variant into a single character.
    Currently available encodings are:
        binary         - reference allele is 0 and any non-reference allele is mutation encoded as 1
        nd16 (default) - full nucleotide diploid model, available in phased and unphased form which
                         is chosen based on the phasing of the VCF file.
    """)
    args = parser.parse_args()
    return args


def write_nexus(names, sequences, file, encoding):
    ntax = len(names)
    nchar = len(sequences[0])
    header = "#nexus\n"
    header += "begin data;\n"
    header += "dimensions ntax=%i nchar=%i;\n" % (ntax, nchar)
    header += "format datatype=\"%s\" gap=- missing=?;\n" % datatype[encoding]
    header += "matrix\n"
    items = [
        "\n".join([name, seq])
        for name, seq in zip(names, sequences)
        ]
    footer = ";\nend;"
    with open(file, "w") as nexus:
        nexus.write(header)
        nexus.write("\n")
        nexus.write("\n\n".join(items))
        nexus.write("\n\n")
        nexus.write(footer)

def vcf2nexus(vcf, nexus, encoding="nd16"):
    names, sequences = parse_vcf(vcf, encoding)
    write_nexus(names, sequences, nexus, encoding)


def parse_vcf(file, encoding):
    with pysam.VariantFile(file) as vcf:
        names = list(vcf.header.samples)
        sequences = get_sequences(vcf.fetch(), encoding)
    return names, sequences


def get_sequences(variants, encoding="nd16"):
    sequences = [variant2bases(variant, encoding) for variant in variants]
    sequences = transpose_list(sequences)
    sequences = map("".join, sequences)
    return list(sequences)


def transpose_list(lst):
    return list(map(list, zip(*lst)))


def variant2bases(variant, encoding="nd16"):
    return [sample2base(sample, encoding) for sample in variant.samples.itervalues()]


def sample2base(sample, encoding="nd16"):
    if encoding == "binary":
        alleles = [
            1 if allele is not None and allele > 1 else allele
            for allele in sample.allele_indices
            ]
    else:
        alleles = sample.alleles
    genome = ["." if allele is None else str(allele) for allele in alleles]
    genome = "".join(genome)
    return translate_genome(genome, encoding, sample.phased)


def translate_genome(genome, encoding="nd16", phased=True):
    if len(genome) != 2:
        raise ValueError(f"Genome must be diploid. The ploidy is: {len(genome)}")
    if encoding == "binary":
        return binary[genome]
    if encoding == "nd16" and phased:
        return nd16phased[genome]
    if encoding == "nd16" and not phased:
        return nd16unphased[genome]
    else:
        raise ValueError(f"Encoding \"{encoding}\" is not supported.")


if __name__ == "__main__":
    args = parse_args()
    vcf2nexus(**vars(args))

# '/usr/bin/env python3
"""Parse a multi-sample VCF file and produce a FASTA file."""
import os
import textwrap
import argparse
import pysam
from Bio import SeqIO

nd16phased = {
    "AA": "0",
    "AC": "1",
    "AG": "2",
    "AT": "3",
    "CA": "4",
    "CC": "5",
    "CG": "6",
    "CT": "7",
    "GA": "8",
    "GC": "9",
    "GG": "a",
    "GT": "b",
    "TA": "c",
    "TC": "d",
    "TG": "e",
    "TT": "f",
    # Phased with one unknown
    ".A": "g",
    ".C": "h",
    ".G": "i",
    ".T": "j",
    "A.": "n",
    "C.": "o",
    "G.": "p",
    "T.": "q",
    "..": "?",
    # haploid states are homozygous genotypes
    "a": "0",
    "A": "0",
    "c": "5",
    "C": "5",
    "g": "a",
    "G": "a",
    "t": "f",
    "T": "f"
}

nd16unphased = {
    "AA": "0",
    "CC": "5",
    "GG": "a",
    "TT": "f",
    "AC": "M",
    "CA": "M",
    "AG": "R",
    "GA": "R",
    "AT": "W",
    "TA": "W",
    "CG": "S",
    "GC": "S",
    "CT": "Y",
    "TC": "Y",
    "GT": "K",
    "TG": "K",
    "..": "?",
    # haploid states are homozygous genotypes
    "a": "0",
    "A": "0",
    "c": "5",
    "C": "5",
    "g": "a",
    "G": "a",
    "t": "f",
    "T": "f"
}

binary = {
    "00": "0",
    "11": "1",
    "01": "1",
    "10": "1",
    ".1": "1",
    "1.": "1",
    "0.": "?",
    ".0": "?",
    "..": "?"
}


def parse_args():
    parser = argparse.ArgumentParser("Transform a multi-sample VCF file into a Fasta format")
    parser.add_argument("vcf", type=str, help="A multi-sample VCF file")
    parser.add_argument("fasta", type=str, help="A path to an output fasta file")
    parser.add_argument("--encoding", type=str, choices=["nd16", "binary"], default="nd16",
                        help="""
    A type of encoding used during to translate a diploid variant into a single character.
    Currently available encodings are:
        binary         - reference allele is 0 and any non-reference allele is mutation encoded as 1
        nd16 (default) - full nucleotide diploid model, available in phased and unphased form which
                         is chosen based on the phasing of the VCF file.
    """)
    parser.add_argument("--ref", type=str,
                        help="A reference Fasta file to fill in genotypes not present in the VCF",
                        default="none")
    parser.add_argument("--debug", type=bool,
                        help="Debug mode (default is false)",
                        default=False)
    args = parser.parse_args()
    return args


def vcf2fasta(vcf, fasta, encoding="nd16", ref="none", debug=False):
    if ref.lower() == "none":
        names, sequences = parse_vcf(vcf, encoding)
        write_fasta(names, sequences, fasta)
    else:
        parse_vcf_ref(vcf, fasta, encoding, ref, debug)


def write_fasta(names, sequences, file):
    items = [
        "\n".join([">" + name, textwrap.fill(seq, width=80)])
        for name, seq in zip(names, sequences)
    ]
    with open(file, "w") as fasta:
        fasta.write("\n\n".join(items))
        fasta.write("\n")


def parse_vcf(file, encoding):
    with pysam.VariantFile(file) as vcf:
        names = list(vcf.header.samples)
        sequences = get_sequences(vcf.fetch(), encoding)
    return names, sequences


def parse_vcf_ref(vcf, fasta, encoding, ref, debug=False):
    if os.path.exists(ref):
        print("using ref fasta: %s" % ref)
        ref_sequence = ""
        # parse ref fasta
        with open(ref) as handle:
            # assumes only one ref sequence in fasta file
            for record in SeqIO.parse(handle, "fasta"):
                ref_sequence_raw = record.seq
        ref_translated = [translate_genome(allele, encoding) for allele in ref_sequence_raw]
        ref_sequence = "".join(ref_translated)
        # parse vcf
        with pysam.VariantFile(vcf) as vcf:
            variants = vcf.fetch()
            samples_map = {}
            for variant in variants:
                for sample in variant.samples.itervalues():
                    base = sample2base(sample, encoding)
                    sample_name = sample.name
                    pos = "%s,%d" % (variant.chrom, variant.pos)
                    if sample_name in samples_map:
                        sample_sequence_map = samples_map.get(sample_name)
                        sample_sequence_map[pos] = base
                    else:
                        sample_sequence_map = {pos: base}
                        samples_map[sample_name] = sample_sequence_map
        if debug:
            print(samples_map)
            print([len(samples_map[i]) for i in samples_map])
        # assume ref sequence is sorted with one-based indexing 1, ..., L
        # open fasta file
        writer = open(fasta, "w")
        for sample_name in samples_map:
            sequence = list(str(ref_sequence))
            variants_map = samples_map[sample_name]
            for variant_pos in variants_map:
                items = variant_pos.split(",")
                chrom = items[0]
                pos = int(items[1]) - 1  # VCF uses one-based index
                sequence[pos] = variants_map[variant_pos]
            sequence_str = "".join(sequence)
            ref_sequence_str = str(ref_sequence)
            if debug:
                variant_index = 555  # one-based position of first variant in 1459.vcf file
                index = variant_index - 1  # zero-based index
                print("ref at pos %d: \n%s" % (index, ref_sequence[index]))
                print("%s at pos %d: \n%s" % (sample_name, index, sequence_str[index]))
                print("ref raw seq: \n%s" % str(ref_sequence_raw)[0:10])
                print("ref: \n%s" % ref_sequence_str[0:10])
                print("%s: \n%s" % (sample_name, sequence_str[0:10]))
            # write sequences line by line
            writer.write(">" + sample_name)
            writer.write("\n")
            writer.write(sequence_str)
            writer.write("\n")
        # close fasta file
        writer.close()
    else:
        # throw exception here
        print("warning: could not find ref fasta: %s" % ref)
    print()


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
    genome = genome.upper()
    # if len(genome) != 2:
    #     raise ValueError(f"Genome must be diploid. The ploidy is: {len(genome)}")
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
    vcf2fasta(**vars(args))
# file = "../data/1459.vcf"
# fasta = file.replace(".vcf", ".fasta")
# encoding = "nd16"
# ref = "../data/filtered_sequence.fna"
# parse_vcf_ref(file, fasta, encoding, ref)

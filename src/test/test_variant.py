import pysam
import pytest
from vcf2fasta import sample2base
from vcf2fasta import variant2bases
from vcf2fasta import get_sequences


@pytest.fixture
def vcf():
    vcf = pysam.VariantHeader()
    vcf.add_sample("sample1")
    vcf.add_sample("sample2")
    vcf.contigs.add("1")
    vcf.contigs.add("2")
    vcf.formats.add("GT", 1, "String", "Genotype")
    return vcf


@pytest.fixture
def variants(vcf):
    return [
        vcf.new_record(
            contig="1",
            alleles=("A", "C"),
            samples=[{"GT" : (0,0), "phased" : True}, {"GT" : (0,1), "phased" : True}]
            ),
        vcf.new_record(
            contig="2",
            alleles=("A", "T"),
            samples=[{"GT" : (0,0), "phased" : True}, {"GT" : (0,1), "phased" : True}]
            )
        ]

@pytest.fixture
def variant(vcf):
    return vcf.new_record(
        contig="1",
        alleles=("A", "C"),
        samples=[{"GT" : (0,0)}, {"GT" : (0,1), "phased" : True}]
        )


@pytest.fixture
def samples(variant):
    return variant.samples.values()


def test_sample2base(variant):
    assert sample2base(variant.samples.values()[0]) == "0"
    assert sample2base(variant.samples.values()[1]) == "1"


def test_variant2bases(variant):
    assert variant2bases(variant) == ["0", "1"]


def test_get_sequences(variants):
    assert get_sequences(variants) == ["00", "13"]

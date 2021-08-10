import pytest
import itertools
from contextlib import contextmanager
from vcf2fasta import translate_genome
#from vcf2fasta import translate_alleles


phased_data_genome = list(map("".join, itertools.product("ACGT", repeat=2))) + [".."]
phased_data_translated = list(map(str, range(0,10))) + list(map(chr, range(97,97+6))) + ["?"]
phased_data = zip(phased_data_genome, phased_data_translated)


@pytest.mark.parametrize("index, value", [(0, "AA"), (1, "AC")])
def test_phased_data_genome(index, value):  
    assert phased_data_genome[index] == value


@pytest.mark.parametrize("index, value", [(0, "0"), (1, "1")])
def test_phased_data_expected(index, value):
    assert phased_data_translated[index] == value


@pytest.mark.parametrize("data, value", [(phased_data_genome, 17), (phased_data_translated, 17)])
def test_phased_data_len(data, value):
    assert len(data) == value


@pytest.mark.parametrize("genome, expected", [
    ("11", "1"),
    ("00", "0"),
    ("01", "1"),
    (".0", "?"),
    ("1.", "1")
    ])
def test_translate_genome_binary(genome, expected):
    assert translate_genome(genome, encoding="binary") == expected


@pytest.mark.parametrize("genome, expected", phased_data)
def test_translate_genome_phased(genome, expected):
    assert translate_genome(genome) == expected


@pytest.mark.parametrize("genome, expected", [
    ("AA", "0"),
    ("CC", "5"),
    ("GG", "a"),
    ("TT", "f"),
    ("..", "?"),
    ])
def test_translate_genome_unphased(genome, expected):
    assert translate_genome(genome, phased=False) == expected


@pytest.mark.parametrize("x", ["A","C","G","T"])
@pytest.mark.parametrize("y", ["A","C","G","T"])
def test_translate_genome_unphased_symetric(x, y):
    translation1 = translate_genome("".join([x,y]), phased=False)
    translation2 = translate_genome("".join([y,x]), phased=False)
    assert translation1 == translation2


@contextmanager
def does_not_raise():
    yield


@pytest.mark.parametrize("genome, expectation", [
    ("A", pytest.raises(ValueError)),
    ("AA", does_not_raise()),
    ("AAA", pytest.raises(ValueError))
    ])
def test_translate_genome_length(genome, expectation):
    with expectation:
        assert translate_genome(genome) is not None


#@pytest.mark.parametrize("alleles, expected", [
#    (("A", "A"), "0"),
#    ((None, None), "?"),
#    ((None, "A"), "g"),
#    (("A", None), "n"),
#    ])
#def test_translate_alleles(alleles, expected):
#    assert translate_alleles(alleles) == expected

from karyopype.karyopype import *
import pandas as pd


def test_get_chromosomes_hg38():
    cs = get_chromsizes("hg38")
    assert len(cs.keys()) == 24


def test_get_chromosomes_nomLee3():
    cs = get_chromsizes("nomLeu3")
    assert len(cs.keys()) == 26


def test_get_chromosomes_X():
    cs = get_chromsizes("mm10")
    assert "chrX" in list(cs.keys())


def test_add_chromsizes_colors():
    cs = get_chromsizes('hg19')
    withcols = add_chromsize_colors(cs, 'green')
    for k in cs:
        assert cs[k][1] == 'green'


def test_parse_regions_None():
    p = parse_regions()
    assert p[0] is None
    assert p[1] is True


def test_parse_regions_str():
    p = parse_regions("data/testRegions/regions.bed")
    assert isinstance(p[0], pd.DataFrame)
    assert p[1] is False


def test_parse_regions_Path():
    p = parse_regions(Path("data/testRegions/regions.bed"))
    assert isinstance(p[0], pd.DataFrame)
    assert p[1] is False

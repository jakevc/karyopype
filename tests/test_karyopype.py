import karyopype.karyopype as kp
import pandas as pd
from pkg_resources import resource_listdir, resource_filename


def test_get_chromosomes_hg38():
    cs = kp.get_chromsizes("hg38")
    assert len(cs.keys()) == 24


def test_get_chromosomes_nomLee3():
    cs = kp.get_chromsizes("nomLeu3")
    assert len(cs.keys()) == 26


def test_get_chromosomes_X():
    cs = kp.get_chromsizes("mm10")
    assert "chrX" in list(cs.keys())


def test_add_chromsizes_colors():
    cs = kp.get_chromsizes('hg19')
    withcols = kp.add_chromsize_colors(cs, 'green')
    for k in cs:
        assert cs[k][1] == 'green'


def test_parse_regions_None():
    p = kp.parse_regions()
    assert p[0] is None
    assert p[1] is True


def test_parse_regions_str():
    dfile = resource_filename(__name__, 'data/testRegions/regions.bed')
    p = kp.parse_regions(dfile)
    assert isinstance(p[0], pd.DataFrame)
    assert p[1] is False


def test_parse_regions_Path():
    dfile = resource_filename(__name__, 'data/testRegions/regions.bed')
    p = kp.parse_regions(dfile)
    assert isinstance(p[0], pd.DataFrame)
    assert p[1] is False


def test_parse_regions_df():
    dat = resource_filename(__name__, 'data/testRegions/regions.bed')
    datdf = pd.read_csv(dat, sep=' ', header=None)
    # test truncation of unncesessary columns
    datdf['extra'] = 'extra'
    p = kp.parse_regions(datdf)
    assert list(p[0].columns) == ['chrom', 'start', 'end']


def test_list_species():
    known = [
            "gorGor3", "hg19", "hg38", "mm10", "nomLeu3", "panTro4", "rheMac8"
            ]
    assert set(kp.list_species()) == set(known)


def test_multiregions():
    r1 = resource_filename(__name__, 'data/testRegions/regions.bed')
    r2 = resource_filename(__name__, 'data/testRegions/regions2.bed')
    a = kp.plot_karyopype("hg19", regions=[r1, r2])
    assert a.__name__  == 'matplotlib.pyplot'

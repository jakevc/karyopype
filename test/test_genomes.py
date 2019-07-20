def test_get_chromosomes_hg38():
    cs = gn.get_chromsizes("hg38")
    assert len(cs.keys()) == 24

def test_get_chromosomes_nomLee3():
    cs = gn.get_chromsizes("nomLeu3")
    assert len(cs.keys()) == 26

def test_get_chromosomes_X():
    cs = gn.get_chromsizes("mm10")
    assert "chrX" in list(cs.keys())

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from overlap_core import *

def f(c, s, e, n="x"): return make_feature(c, s, e, n)

def test_overlap_basic():
    assert features_overlap(f("chr1", 0, 10), f("chr1", 5, 15))
    assert overlap_length(f("chr1", 0, 10), f("chr1", 5, 15)) == 5

def test_touching_is_not_overlap():          # half-open BED coordinates
    assert not features_overlap(f("chr1", 0, 10), f("chr1", 10, 20))

def test_different_chrom():
    assert not features_overlap(f("chr1", 0, 10), f("chr2", 0, 10))

def test_snp_in_gene():
    assert overlap_length(f("chr1", 1000, 5000), f("chr1", 2000, 2001)) == 1

def test_parse_and_errors():
    feats, errs = parse_bed_text("chr1\t1\t5\tA\nbad line\nchr1\t9\t3\tB\n#c\n")
    assert len(feats) == 1 and len(errs) == 2

def test_sample_data():
    g, _ = read_feature_file("data/genes.bed"); s, _ = read_feature_file("data/snps.bed")
    ov = find_all_overlaps(g, s)
    assert len(ov) == 4
    assert features_without_overlap(g, ov) == []

def test_sequence_info():
    g, _ = read_feature_file("data/genes.bed"); s, _ = read_feature_file("data/snps.bed")
    df = overlaps_to_dataframe(find_all_overlaps(g, s))
    with open("data/reference.fa") as h: ref = load_reference_fasta(h)
    out = add_sequence_info(df, ref)
    assert "GC_percent" in out.columns

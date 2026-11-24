"""
Two-enzyme cloning: the two ends of the vector are cut by different enzymes, possibly with
3' overhangs, palindromic or degenerate sites. Tested on a real vector (pMCB306, BstXI + BlpI)
and, independently of the code's own cut arithmetic, on synthetic vectors for many enzyme pairs.
"""

import itertools
import random
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import sgrna_oligo_designer as sod

FIXTURES = Path(__file__).resolve().parent.parent / "PyShinyApp" / "vectors"
SG = "GACCAGGATGGGCACCACCC"          # the placeholder guide that ships in pMCB306


@pytest.fixture(scope="module")
def pmcb306():
    return sod.read_vector_gb(FIXTURES / "pMCB306.gb")


# --- real vector: pMCB306 (BstXI 3' overhang on the left, BlpI 5' overhang on the right) -------

def test_pmcb306_geometry_and_flanks(pmcb306):
    site = sod.find_cloning_site(pmcb306, "BstXI", enzyme_b="BlpI")
    assert site.orientation == "forward" and site.enzymes == "BstXI + BlpI"
    assert (site.left_overhang, site.right_overhang) == ("CTTG", "TAA")       # 4-nt 3' / 3-nt 5'
    assert (site.left_extra, site.right_extra) == ("TTG", "GTTTAAGAGC")
    assert not site.warnings


def test_pmcb306_oligos_follow_the_bstxi_blpi_layout(pmcb306):
    site = sod.find_cloning_site(pmcb306, "BstXI", enzyme_b="BlpI")
    (row,) = sod.add_cloning_oligos([{"sgRNA_id": "a", "sgRNA_seq": SG}], site, g_mode="if_needed")
    assert row["oligo_sense"] == "TTG" + SG + "GTTTAAGAGC"
    assert row["oligo_antisense"] == "TTAGCTCTTAAAC" + sod.reverse_complement(SG) + "CAACAAG"


def test_pmcb306_cloning_its_own_placeholder_restores_the_vector(pmcb306):
    site = sod.find_cloning_site(pmcb306, "BstXI", enzyme_b="BlpI")
    rows = sod.add_cloning_oligos([{"sgRNA_id": "nt", "sgRNA_seq": SG}], site, g_mode="never")
    (out,) = sod.build_plasmid_records(rows, pmcb306, site)
    assert str(out.seq).upper() == str(pmcb306.seq).upper()


def test_pmcb306_enzyme_order_does_not_matter_and_both_are_single_cutters(pmcb306):
    a = sod.find_cloning_site(pmcb306, "BstXI", enzyme_b="BlpI")
    b = sod.find_cloning_site(pmcb306, "BlpI", enzyme_b="BstXI")
    assert (a.left_overhang, a.right_overhang, a.left_extra) == (b.left_overhang, b.right_overhang, b.left_extra)
    assert {"BstXI", "BlpI"} <= set(sod.single_cutters(pmcb306))


def test_pmcb306_spacer_check_ignores_the_regenerated_vector_sites(pmcb306):
    site = sod.find_cloning_site(pmcb306, "BstXI", enzyme_b="BlpI")
    rows = sod.add_cloning_oligos([
        {"sgRNA_id": "ok", "sgRNA_seq": SG},
        {"sgRNA_id": "bad", "sgRNA_seq": "GCTAAGCAAAAGGCACCACC"},     # contains a BlpI site GCTNAGC
    ], site, g_mode="if_needed")
    assert rows[0]["warnings"] == []
    assert any("BlpI" in w and "recognition site" in w for w in rows[1]["warnings"])


def test_user_flanks_override_the_suggestion(pmcb306):
    site = sod.find_cloning_site(pmcb306, "BstXI", enzyme_b="BlpI", left_extra="tt g", right_extra="")
    assert (site.left_extra, site.right_extra) == ("TTG", "")
    with pytest.raises(ValueError, match="only A, C, G and T"):
        sod.find_cloning_site(pmcb306, "BstXI", enzyme_b="BlpI", left_extra="TTX")
    with pytest.raises(ValueError, match="two-enzyme mode"):
        sod.find_cloning_site(pmcb306, "BsmBI", left_extra="TTG")


def test_two_enzyme_error_messages(pmcb306):
    seq = str(pmcb306.seq).upper()
    cand = sod.two_enzyme_candidates()
    never = next(n for n, e in cand.items() if len(sod._cuts(seq, e)) == 0)
    with pytest.raises(ValueError, match=f"{never} cuts this vector 0 time"):
        sod.find_cloning_site(pmcb306, never, enzyme_b="BlpI")
    many = next(n for n, e in cand.items() if len(sod._cuts(seq, e)) > 1)
    with pytest.raises(ValueError, match=f"{many} cuts this vector [2-9]\\d* time"):
        sod.find_cloning_site(pmcb306, "BstXI", enzyme_b=many)
    with pytest.raises(ValueError, match="two different enzymes"):
        sod.find_cloning_site(pmcb306, "BlpI", enzyme_b="BlpI")
    with pytest.raises(ValueError, match="Unsupported enzyme"):
        sod.find_cloning_site(pmcb306, "BlpI", enzyme_b="NotAnEnzyme")


# --- synthetic vectors: independent check across many enzyme pairs ---------------------------------

def _instance(site: str) -> str:
    """A concrete sequence matching an IUPAC site (N and friends -> fixed bases)."""
    fill = {"N": "A", "R": "A", "Y": "C", "S": "C", "W": "A", "K": "G", "M": "A", "B": "C",
            "D": "A", "H": "A", "V": "A"}
    return "".join(fill.get(c, c) for c in site)


def _filler(n, seed, avoid):
    rng = random.Random(seed)
    while True:
        s = "".join(rng.choice("ACGT") for _ in range(n))
        if not any(p in s for p in avoid):
            return s


PAIRS = [("EcoRI", "AgeI"), ("BamHI", "EcoRI"), ("XhoI", "EcoRI"), ("HindIII", "BglII"),
         ("BstXI", "BlpI"), ("BsmBI", "BsaI"), ("SapI", "AgeI"), ("KpnI", "SacI"),
         ("PstI", "EcoRI"), ("AgeI", "BlpI")]


def _synthetic(name_a, name_b):
    cand = sod.two_enzyme_candidates()
    ea, eb = cand[name_a], cand[name_b]
    sa, sb = _instance(str(ea.site)), _instance(str(eb.site))
    avoid = {sa, sb, sod.reverse_complement(sa), sod.reverse_complement(sb)}
    for seed in range(80):
        seq = (_filler(250, 3 * seed + 1, avoid) + sa + _filler(40, 3 * seed + 2, avoid) + sb
               + _filler(250, 3 * seed + 3, avoid))
        if len(sod._cuts(seq, ea)) == 1 and len(sod._cuts(seq, eb)) == 1:
            return SeqRecord(Seq(seq), id="V", name="V", annotations={"molecule_type": "DNA", "topology": "circular"})
    raise AssertionError(f"no clean synthetic vector for {name_a}/{name_b}")


@pytest.mark.parametrize("name_a,name_b", PAIRS)
def test_cloning_back_the_original_stuffer_releases_exactly_the_digest_strands(name_a, name_b):
    rec = _synthetic(name_a, name_b)
    seq = str(rec.seq)
    site = sod.find_cloning_site(rec, name_a, "forward", enzyme_b=name_b, left_extra="", right_extra="")

    # what the digest releases, from Biopython's own cut positions (independent of the trim arithmetic)
    cand = sod.two_enzyme_candidates()
    (ta, ba), = sod._cuts(seq, cand[name_a])
    (tb, bb), = sod._cuts(seq, cand[name_b])
    (tl, bl), (tr, br) = sorted([(ta, ba), (tb, bb)], key=lambda c: min(c))
    released_top, released_bottom = seq[tl:tr], sod.reverse_complement(seq[bl:br])

    original_core = seq[site.replace_start:site.replace_end]
    top, bottom = sod.design_oligos(original_core, site)
    assert (top, bottom) == (released_top, released_bottom)

    # and putting that core back reproduces the vector exactly
    (out,) = [sod.insert_sgrna(rec, site, original_core, original_core, "x")]
    assert str(out.seq) == seq


@pytest.mark.parametrize("name_a,name_b", PAIRS)
def test_new_guide_gives_a_plasmid_that_contains_both_oligos_in_place(name_a, name_b):
    rec = _synthetic(name_a, name_b)
    site = sod.find_cloning_site(rec, name_a, "forward", enzyme_b=name_b, left_extra="", right_extra="")
    rows = sod.add_cloning_oligos([{"sgRNA_id": "a", "sgRNA_seq": "ACGTTGCAAGCTTGGCATCA"}], site, g_mode="never")
    (out,) = sod.build_plasmid_records(rows, rec, site)
    plasmid = str(out.seq)
    top, bottom = rows[0]["oligo_sense"], rows[0]["oligo_antisense"]
    assert top in plasmid                                           # top oligo sits in the top strand
    assert sod.reverse_complement(bottom) in plasmid                # bottom oligo pairs with it
    assert len(plasmid) == len(rec.seq) - site.removed_length + len(rows[0]["insert_core"])


def test_identical_overhangs_warn():
    """BamHI and BglII leave the same GATC overhang, so the insert could go in either way round."""
    rec = _synthetic("BamHI", "BglII")
    site = sod.find_cloning_site(rec, "BamHI", "forward", enzyme_b="BglII")
    assert any("identical" in w for w in site.warnings)


def test_two_enzyme_candidates_include_the_classics():
    cand = sod.two_enzyme_candidates()
    assert {"EcoRI", "AgeI", "BamHI", "XhoI", "HindIII", "BglII", "BstXI", "BlpI", "BsmBI"} <= set(cand)
    assert "SmaI" not in cand and "EcoRV" not in cand          # blunt cutters: no overhang

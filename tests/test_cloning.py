"""
Tests for the enzyme-driven design: in-silico digestion, overhang derivation,
both site layouts, both orientations, and the error cases. They use synthetic
vectors so every enzyme can be tested (the example vectors are BsmBI-only).
"""

import random

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

import sgrna_oligo_designer as sod

OL, OR = "CACC", "GTTT"
SG = "GATTCTAAAACGGATTACCA"
ENZ_PARAMS = ["BsmBI", "Esp3I", "BsaI", "BbsI"]


def _filler(n, seed):
    """Random ACGT with no site for any supported enzyme, on either strand."""
    pats = {p for e in sod.ENZYMES
            for p in (sod.get_enzyme(e).site, sod.reverse_complement(sod.get_enzyme(e).site))}
    rng = random.Random(seed)
    while True:
        s = "".join(rng.choice("ACGT") for _ in range(n))
        if not any(p in s for p in pats):
            return s


BB_L, BB_R, MID = _filler(300, 1), _filler(300, 2), _filler(40, 3)


def make_vector(enzyme_name, layout, left_ov=OL, right_ov=OR, mid=MID, bb_l=BB_L, bb_r=BB_R,
                g_in_vector=False, topology="circular"):
    """
    Build a circular vector with two sites around a stuffer.

    layout 'A': sites sit inside the removed piece (lost on cutting, lentiCRISPR-like)
    layout 'B': sites stay in the backbone, facing the stuffer
    Returns (record, expected_plasmid_core_context) where the context is
    (left_part, right_part) of the final sequence around the insert.
    """
    enz = sod.get_enzyme(enzyme_name)
    S, R = enz.site, sod.reverse_complement(enz.site)
    gap = "A" * (enz.cut_offset - len(S))          # the N bases between site and cut
    if g_in_vector:                                 # like ipUSEPR: CACC|G GAGACG... the G is that N
        gap = "G" + gap[1:]
    if layout == "A":
        seq = bb_l + left_ov + gap + R + "T" + mid + "T" + S + gap + right_ov + bb_r
        ctx = (bb_l + left_ov, right_ov + bb_r)
    else:
        seq = bb_l + "AT" + S + gap + left_ov + mid + right_ov + gap + R + "AT" + bb_r
        ctx = (bb_l + "AT" + S + gap + left_ov, right_ov + gap + R + "AT" + bb_r)
    rec = SeqRecord(Seq(seq), id="V", name="V", description="synthetic",
                    annotations={"molecule_type": "DNA", "topology": topology})
    return rec, ctx


def rc_record(rec):
    return rec.reverse_complement(**sod._RC_KEEP)


# --- digestion: overhangs, both layouts, all enzymes -----------------------------

@pytest.mark.parametrize("layout", ["A", "B"])
@pytest.mark.parametrize("enzyme", ENZ_PARAMS)
def test_overhangs_read_from_the_cut_not_typed(enzyme, layout):
    rec, _ = make_vector(enzyme, layout, left_ov="ACCG", right_ov="GTTT")
    site = sod.find_cloning_site(rec, enzyme)
    assert (site.left_overhang, site.right_overhang) == ("ACCG", "GTTT")
    assert site.orientation == "forward"


@pytest.mark.parametrize("layout", ["A", "B"])
@pytest.mark.parametrize("enzyme", ENZ_PARAMS)
def test_plasmid_is_backbone_plus_oligo_duplex(enzyme, layout):
    rec, (left, right) = make_vector(enzyme, layout)
    site = sod.find_cloning_site(rec, enzyme)
    rows = sod.add_cloning_oligos([{"sgRNA_id": "a", "sgRNA_seq": SG}], site)
    (out,) = sod.build_plasmid_records(rows, rec, site)
    assert str(out.seq) == left + "G" + SG + right
    # oligos anneal into exactly this: sense = OL+G+sg ; antisense = rc(G+sg+OR)
    assert rows[0]["oligo_sense"] == OL + "G" + SG
    assert rows[0]["oligo_antisense"] == sod.reverse_complement("G" + SG + OR)


def test_layout_a_loses_sites_layout_b_keeps_them():
    for layout, remaining in (("A", 0), ("B", 2)):
        rec, _ = make_vector("BsaI", layout)
        site = sod.find_cloning_site(rec, "BsaI")
        (out,) = sod.build_plasmid_records(
            sod.add_cloning_oligos([{"sgRNA_id": "a", "sgRNA_seq": SG}], site), rec, site)
        assert len(sod._find_sites(str(out.seq), sod.get_enzyme("BsaI"))) == remaining


def test_reverse_orientation_is_detected_and_round_trips():
    fwd_rec, _ = make_vector("BsaI", "A")
    rev_rec = rc_record(fwd_rec)
    site = sod.find_cloning_site(rev_rec, "BsaI")
    assert site.orientation == "reverse"
    assert (site.left_overhang, site.right_overhang) == (OL, OR)

    rows = sod.add_cloning_oligos([{"sgRNA_id": "a", "sgRNA_seq": SG}], site)
    (rev_out,) = sod.build_plasmid_records(rows, rev_rec, site)
    fsite = sod.find_cloning_site(fwd_rec, "BsaI")
    frows = sod.add_cloning_oligos([{"sgRNA_id": "a", "sgRNA_seq": SG}], fsite)
    (fwd_out,) = sod.build_plasmid_records(frows, fwd_rec, fsite)
    assert str(rev_out.seq) == sod.reverse_complement(str(fwd_out.seq))
    assert [r["oligo_sense"] for r in rows] == [r["oligo_sense"] for r in frows]


def test_forced_orientation_overrides_auto():
    rec, _ = make_vector("BsaI", "A", left_ov="TTTT", right_ov="CCCC")      # no scaffold overhang
    auto = sod.find_cloning_site(rec, "BsaI")
    assert auto.orientation == "forward" and "Could not confirm" in auto.warnings[0]
    assert sod.find_cloning_site(rec, "BsaI", "forward").warnings == ()
    rev = sod.find_cloning_site(rec, "BsaI", "reverse")
    assert (rev.left_overhang, rev.right_overhang) == ("GGGG", "AAAA")


# --- errors ---------------------------------------------------------------------

def test_wrong_enzyme_suggests_the_right_one():
    rec = sod.read_vector_gb(
        __import__("pathlib").Path(__file__).resolve().parent.parent / "PyShinyApp" / "vectors" / "ipUSEPR.gb")
    with pytest.raises(ValueError) as e:
        sod.find_cloning_site(rec, "BsaI")
    msg = str(e.value)
    assert "BsaI" in msg and "6 time(s)" in msg and "BsmBI" in msg


def test_enzyme_absent_from_vector():
    rec, _ = make_vector("BsaI", "A")
    with pytest.raises(ValueError, match="BbsI .* 0 time"):
        sod.find_cloning_site(rec, "BbsI")


def test_unknown_enzyme_and_orientation():
    rec, _ = make_vector("BsaI", "A")
    with pytest.raises(ValueError, match="Unsupported enzyme"):
        sod.find_cloning_site(rec, "EcoRI")
    with pytest.raises(ValueError, match="orientation"):
        sod.find_cloning_site(rec, "BsaI", "sideways")


def test_linear_vector_rejected():
    rec, _ = make_vector("BsaI", "A", topology="linear")
    with pytest.raises(ValueError, match="circular"):
        sod.find_cloning_site(rec, "BsaI")


def test_cloning_site_spanning_the_origin_is_rejected():
    rec, _ = make_vector("BsaI", "A")
    cut = len(BB_L) + 20                              # inside the removed piece
    rotated = SeqRecord(rec.seq[cut:] + rec.seq[:cut], id="V", name="V",
                        annotations=dict(rec.annotations))
    with pytest.raises(ValueError, match="origin"):
        sod.find_cloning_site(rotated, "BsaI")


def test_site_for_a_different_vector_is_rejected():
    rec, _ = make_vector("BsaI", "A")
    site = sod.find_cloning_site(rec, "BsaI")
    other, _ = make_vector("BsaI", "A", bb_l=BB_L + "ACGT")
    with pytest.raises(ValueError, match="different vector"):
        sod.insert_sgrna(other, site, "G" + SG, SG, "x")


# --- the 5' G -----------------------------------------------------------------------

@pytest.mark.parametrize("sg,mode,core", [
    (SG, "always", "G" + SG), (SG, "if_needed", SG), (SG, "never", SG),
    ("CATTCTAAAACGGATTACCA", "always", "GCATTCTAAAACGGATTACCA"),
    ("CATTCTAAAACGGATTACCA", "if_needed", "GCATTCTAAAACGGATTACCA"),
    ("CATTCTAAAACGGATTACCA", "never", "CATTCTAAAACGGATTACCA"),
])
def test_g_modes(sg, mode, core):
    assert sod.build_insert_core(sg, mode) == core


def test_invalid_g_mode():
    with pytest.raises(ValueError):
        sod.build_insert_core(SG, "sometimes")


@pytest.mark.parametrize("mode", sod.G_MODES)
def test_plasmid_matches_oligos_for_every_g_mode(mode):
    rec, (left, right) = make_vector("BsmBI", "A", g_in_vector=True)
    site = sod.find_cloning_site(rec, "BsmBI")
    rows = sod.add_cloning_oligos([{"sgRNA_id": "a", "sgRNA_seq": SG}], site, g_mode=mode)
    (out,) = sod.build_plasmid_records(rows, rec, site)
    core = rows[0]["insert_core"]
    assert str(out.seq) == left + core + right
    f = next(f for f in out.features if f.qualifiers.get("label") == ["a"])
    assert str(f.extract(out.seq)) == SG


# --- warnings ---------------------------------------------------------------------

@pytest.mark.parametrize("sg", [
    "GG" + "CGTCTC" + "A" * 12,                       # site in the spacer
    "GG" + "GAGACG" + "A" * 12,                       # site on the other strand
    "T" * 15 + "GAGAC",                               # site created at the junction with GTTT
])
def test_spacer_containing_the_enzyme_site_warns(sg):
    rec, _ = make_vector("BsmBI", "A")
    site = sod.find_cloning_site(rec, "BsmBI")
    (row,) = sod.add_cloning_oligos([{"sgRNA_id": "bad", "sgRNA_seq": sg}], site)
    assert any("BsmBI recognition site" in w for w in row["warnings"])


def test_clean_spacer_has_no_warnings_and_wrong_length_warns():
    rec, _ = make_vector("BsmBI", "A")
    site = sod.find_cloning_site(rec, "BsmBI")
    rows = sod.add_cloning_oligos([{"sgRNA_id": "ok", "sgRNA_seq": SG},
                                   {"sgRNA_id": "short", "sgRNA_seq": "ACGTACGTAC"},
                                   {"sgRNA_id": "pam", "sgRNA_seq": SG + "TGG"}], site)
    assert rows[0]["warnings"] == []
    assert "is 10 nt" in rows[1]["warnings"][0]
    assert "NGG PAM" in rows[2]["warnings"][-1]


# --- feature handling ---------------------------------------------------------------

def _feat(start, end, label, strand=1):
    return SeqFeature(FeatureLocation(start, end, strand), type="misc_feature",
                      qualifiers={"label": [label]})


def test_features_inside_removed_dropped_overlapping_clipped_downstream_shifted():
    rec, _ = make_vector("BsaI", "A")
    site = sod.find_cloning_site(rec, "BsaI")
    r0, r1 = site.replace_start, site.replace_end
    rec.features = [
        _feat(10, 30, "upstream"),
        _feat(r0 + 3, r1 - 3, "inside"),
        _feat(0, len(rec.seq), "source_like"),
        _feat(r0 - 10, r0 + 5, "clip_left"),
        _feat(r1 - 5, r1 + 10, "clip_right"),
        _feat(r1 + 20, r1 + 40, "downstream", -1),
    ]
    rows = sod.add_cloning_oligos([{"sgRNA_id": "a", "sgRNA_seq": SG}], site)
    (out,) = sod.build_plasmid_records(rows, rec, site)
    by = {f.qualifiers["label"][0]: f for f in out.features}
    n_new = len(out.seq)
    assert "inside" not in by
    assert (int(by["upstream"].location.start), int(by["upstream"].location.end)) == (10, 30)
    assert (int(by["source_like"].location.start), int(by["source_like"].location.end)) == (0, n_new)
    assert int(by["clip_left"].location.end) == r0
    assert int(by["clip_right"].location.start) == r0 + len(rows[0]["insert_core"])
    d = by["downstream"]
    shift = n_new - len(rec.seq)
    assert (int(d.location.start), int(d.location.end), d.location.strand) == (r1 + 20 + shift, r1 + 40 + shift, -1)


def test_compound_feature_is_remapped_part_by_part():
    rec, _ = make_vector("BsaI", "A")
    site = sod.find_cloning_site(rec, "BsaI")
    r1 = site.replace_end
    rec.features = [SeqFeature(CompoundLocation([FeatureLocation(5, 25, 1),
                                                 FeatureLocation(r1 + 10, r1 + 30, 1)]),
                               type="CDS", qualifiers={"label": ["joined"]})]
    rows = sod.add_cloning_oligos([{"sgRNA_id": "a", "sgRNA_seq": SG}], site)
    (out,) = sod.build_plasmid_records(rows, rec, site)
    f = next(f for f in out.features if f.type == "CDS")
    shift = len(out.seq) - len(rec.seq)
    assert [(int(p.start), int(p.end)) for p in f.location.parts] == [(5, 25), (r1 + 10 + shift, r1 + 30 + shift)]


def test_vector_provided_g_keeps_its_annotation():
    rec, _ = make_vector("BsmBI", "A", g_in_vector=True)
    site = sod.find_cloning_site(rec, "BsmBI")
    g_pos = site.replace_start
    assert str(rec.seq)[g_pos] == "G"
    rec.features = [_feat(g_pos, g_pos + 1, "tss")]
    rows = sod.add_cloning_oligos([{"sgRNA_id": "a", "sgRNA_seq": SG}], site)       # always -> G
    (out,) = sod.build_plasmid_records(rows, rec, site)
    tss = next(f for f in out.features if f.qualifiers["label"] == ["tss"])
    assert (int(tss.location.start), int(tss.location.end)) == (g_pos, g_pos + 1)


def test_input_vector_is_not_modified():
    rec, _ = make_vector("BsaI", "A")
    before = (str(rec.seq), len(rec.features))
    site = sod.find_cloning_site(rec, "BsaI")
    sod.build_plasmid_records(
        sod.add_cloning_oligos([{"sgRNA_id": "a", "sgRNA_seq": SG}], site), rec, site)
    assert (str(rec.seq), len(rec.features)) == before


# --- review follow-ups: warnings, messages, qualifiers ---------------------------

def test_tttt_terminator_in_spacer_warns():
    rec, _ = make_vector("BsmBI", "A")
    site = sod.find_cloning_site(rec, "BsmBI")
    rows = sod.add_cloning_oligos([{"sgRNA_id": "t", "sgRNA_seq": "AAAATTTTAAAACCCCGGGG"},
                                   {"sgRNA_id": "ok", "sgRNA_seq": "AAAATTTAAAACCCCGGGGA"}], site)
    assert any("TTTT" in w and "RNA polymerase III" in w for w in rows[0]["warnings"])
    assert not any("TTTT" in w for w in rows[1]["warnings"])


def test_identical_and_palindromic_overhangs_warn():
    same, _ = make_vector("BsmBI", "A", left_ov="CACC", right_ov="CACC")
    assert any("identical" in w for w in sod.find_cloning_site(same, "BsmBI", "forward").warnings)
    pal, _ = make_vector("BsmBI", "A", left_ov="AATT", right_ov="GTTT")
    assert any("palindromic" in w for w in sod.find_cloning_site(pal, "BsmBI").warnings)
    ok, _ = make_vector("BsmBI", "A")
    assert sod.find_cloning_site(ok, "BsmBI").warnings == ()


def test_large_stuffer_error_reports_both_piece_sizes():
    rec, _ = make_vector("BsmBI", "A", mid=_filler(800, 7), bb_l=_filler(100, 8), bb_r=_filler(100, 9))
    with pytest.raises(ValueError, match=r"\d+ bp piece between the sites and a \d+ bp piece that crosses the origin"):
        sod.find_cloning_site(rec, "BsmBI")


def test_translation_qualifier_dropped_only_when_the_feature_is_edited():
    rec, _ = make_vector("BsaI", "A")
    site = sod.find_cloning_site(rec, "BsaI")
    r0, r1 = site.replace_start, site.replace_end
    cds = lambda a, b, label: SeqFeature(FeatureLocation(a, b, 1), type="CDS",
                                         qualifiers={"label": [label], "translation": ["MKV"]})
    rec.features = [cds(10, 40, "before"), cds(r0 - 20, r1 + 20, "spanning"), cds(r1 + 30, r1 + 60, "after")]
    rows = sod.add_cloning_oligos([{"sgRNA_id": "a", "sgRNA_seq": SG}], site)
    (out,) = sod.build_plasmid_records(rows, rec, site)
    by = {f.qualifiers["label"][0]: f for f in out.features if f.type == "CDS"}
    assert "translation" in by["before"].qualifiers and "translation" in by["after"].qualifiers
    assert "translation" not in by["spanning"].qualifiers


def test_read_vector_gb_gives_friendly_errors(tmp_path):
    from pathlib import Path
    good = (Path(__file__).resolve().parent.parent / "PyShinyApp" / "vectors" / "ipUSEPR.gb").read_text()
    two = tmp_path / "two.gb"
    two.write_text(good + good)
    with pytest.raises(ValueError, match="more than one record"):
        sod.read_vector_gb(two)
    junk = tmp_path / "junk.gb"
    junk.write_text("this is not a genbank file\n")
    with pytest.raises(ValueError, match="Could not read this file as a GenBank"):
        sod.read_vector_gb(junk)


# --- any eligible enzyme ("Other…") ------------------------------------------------------

def _filler_avoiding(n, seed, avoid):
    rng = random.Random(seed)
    while True:
        s = "".join(rng.choice("ACGT") for _ in range(n))
        if not any(p in s for p in avoid):
            return s


def _vector_for(name, layout):
    """Synthetic vector for any eligible enzyme, with overhangs of that enzyme's length."""
    enz = sod.get_enzyme(name)
    k = enz.overhang_len
    left, right = "CACCG"[:k], "GTTTA"[:k]
    avoid = (enz.site, sod.reverse_complement(enz.site))
    for seed in range(60):
        rec, ctx = make_vector(
            name, layout, left_ov=left, right_ov=right,
            bb_l=_filler_avoiding(200, 3 * seed + 1, avoid), bb_r=_filler_avoiding(200, 3 * seed + 2, avoid),
            mid=_filler_avoiding(30, 3 * seed + 3, avoid))
        if len(sod._find_sites(str(rec.seq).upper(), enz)) == 2:
            return rec, ctx, left, right
    raise AssertionError(f"could not build a clean synthetic vector for {name}")


@pytest.mark.parametrize("layout", ["A", "B"])
@pytest.mark.parametrize("name", sorted(sod.all_enzymes()))
def test_every_eligible_enzyme_round_trips(name, layout):
    rec, (left_ctx, right_ctx), left, right = _vector_for(name, layout)
    site = sod.find_cloning_site(rec, name)
    assert (site.left_overhang, site.right_overhang, site.orientation) == (left, right, "forward")
    rows = sod.add_cloning_oligos([{"sgRNA_id": "a", "sgRNA_seq": SG}], site)
    (out,) = sod.build_plasmid_records(rows, rec, site)
    assert str(out.seq) == left_ctx + "G" + SG + right_ctx
    assert rows[0]["oligo_sense"] == left + "G" + SG
    assert rows[0]["oligo_antisense"] == sod.reverse_complement("G" + SG + right)


def test_eligible_list_contents():
    names = set(sod.all_enzymes())
    assert {"BsmBI", "Esp3I", "BsaI", "BbsI", "SapI", "BspQI", "LguI", "BpiI", "Eco31I"} <= names
    assert sod.get_enzyme("SapI").overhang_len == 3
    assert not names & {"EcoRI", "FspEI", "AlwI", "PleI", "BcgI", "BsaXI"}      # palindromic / 1-2 nt / two-sided
    assert all(sod.MIN_OVERHANG <= e.overhang_len <= sod.MAX_OVERHANG for e in sod.all_enzymes().values())
    for bad in ("EcoRI", "AlwI", "NotAnEnzyme"):
        with pytest.raises(ValueError, match="Unsupported enzyme"):
            sod.get_enzyme(bad)


def test_three_nt_overhang_orientation_is_detected_both_ways():
    fwd, _, _, _ = _vector_for("SapI", "A")
    assert sod.find_cloning_site(fwd, "SapI").orientation == "forward"
    assert sod.find_cloning_site(fwd, "SapI").right_overhang == "GTT"
    rev = rc_record(fwd)
    site = sod.find_cloning_site(rev, "SapI")
    assert (site.orientation, site.left_overhang, site.right_overhang) == ("reverse", "CAC", "GTT")


def test_wrong_enzyme_hint_lists_other_enzymes_that_do_cut_twice():
    rec, _, _, _ = _vector_for("SapI", "A")
    with pytest.raises(ValueError) as e:
        sod.find_cloning_site(rec, "BsmBI")
    assert "SapI" in str(e.value) and "exactly two sites for" in str(e.value)

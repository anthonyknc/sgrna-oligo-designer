"""Unit tests for core logic using small synthetic inputs."""

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

import sgrna_oligo_designer as sod

SG = "GATTCTAAAACGGATTACCA"
STUFFER = "GAGACGTGACCGTCTCT"          # 17 nt, non-palindromic
LEFT, RIGHT = "ACTGACTGAAGCTTGCATGC", "TTGGCCAATTCCGGAATTCA"


def _vector(features=()):
    rec = SeqRecord(Seq(LEFT + STUFFER + RIGHT), id="V", name="V",
                    annotations={"molecule_type": "DNA", "topology": "circular"})
    rec.features = list(features)
    return rec


# --- oligo design ---------------------------------------------------------

def test_reverse_complement():
    assert sod.reverse_complement("AACG") == "CGTT"


def test_default_oligos():
    (row,) = sod.add_cloning_oligos([{"sgRNA_id": "a", "sgRNA_seq": SG}])
    assert row["oligo_sense"] == "CACCG" + SG
    assert row["oligo_antisense"] == "AAAC" + sod.reverse_complement(SG) + "C"


def test_custom_prefixes():
    (row,) = sod.add_cloning_oligos([{"sgRNA_id": "a", "sgRNA_seq": SG}],
                                    top_prefix="CACC", bot_prefix="AAAC", bot_suffix="")
    assert row["oligo_sense"] == "CACC" + SG
    assert row["oligo_antisense"] == "AAAC" + sod.reverse_complement(SG)


def test_whitespace_and_case_cleaned():
    messy = " gattctaaaa\n cggattacca\t"
    (row,) = sod.add_cloning_oligos([{"sgRNA_id": "a", "sgRNA_seq": messy}])
    assert row["sgRNA_seq"] == SG


@pytest.mark.parametrize("bad", ["", "   ", "GATTCTAAAACGGATTACCN", "GATTCTAAAACGGATTAC-A"])
def test_invalid_sgrna_raises(bad):
    with pytest.raises(ValueError):
        sod.add_cloning_oligos([{"sgRNA_id": "a", "sgRNA_seq": bad}])


def test_wrong_length_and_pam_warn(capsys):
    sod.add_cloning_oligos([{"sgRNA_id": "short", "sgRNA_seq": "ACGTACGTAC"}])
    sod.add_cloning_oligos([{"sgRNA_id": "pam", "sgRNA_seq": SG + "TGG"}])
    err = capsys.readouterr().err
    assert "'short' is 10 nt" in err
    assert "NGG PAM" in err


def test_missing_column_raises():
    with pytest.raises(ValueError, match="not found"):
        sod.add_cloning_oligos([{"id": "a", "seq": SG}])


def test_idt_rows_interleaved():
    rows = sod.add_cloning_oligos([{"sgRNA_id": "a", "sgRNA_seq": SG},
                                   {"sgRNA_id": "b", "sgRNA_seq": SG}])
    idt = sod.build_idt_rows(rows)
    assert [r["Name"] for r in idt] == ["a_s", "a_a", "b_s", "b_a"]
    assert {r["Scale"] for r in idt} == {"25nm"} and {r["Purification"] for r in idt} == {"STD"}


# --- stuffer location -----------------------------------------------------

def test_stuffer_forward():
    assert sod.validate_stuffer(_vector(), STUFFER) == (len(LEFT), +1)


def test_stuffer_reverse_strand():
    rec = SeqRecord(Seq(LEFT + sod.reverse_complement(STUFFER) + RIGHT), id="V", name="V")
    assert sod.validate_stuffer(rec, STUFFER) == (len(LEFT), -1)


def test_stuffer_errors():
    with pytest.raises(ValueError, match="not found"):
        sod.validate_stuffer(_vector(), "CCCCCCCCCCCCCCCCCCCA")
    with pytest.raises(ValueError, match="palindromic"):
        sod.validate_stuffer(_vector(), "GAATTC")
    with pytest.raises(ValueError, match="empty"):
        sod.validate_stuffer(_vector(), "")
    with pytest.raises(ValueError, match="non-ACGT"):
        sod.validate_stuffer(_vector(), "GAGACGN")
    dup = SeqRecord(Seq(LEFT + STUFFER + RIGHT + STUFFER), id="V", name="V")
    with pytest.raises(ValueError, match="2 times"):
        sod.validate_stuffer(dup, STUFFER)


# --- insertion ------------------------------------------------------------

def test_insert_forward_and_reverse_strand():
    fwd = sod.insert_sgrna(_vector(), STUFFER, SG, "sgA")
    assert str(fwd.seq) == LEFT + SG + RIGHT
    rev_vec = SeqRecord(Seq(LEFT + sod.reverse_complement(STUFFER) + RIGHT), id="V", name="V",
                        annotations={"molecule_type": "DNA"})
    rev = sod.insert_sgrna(rev_vec, STUFFER, SG, "sgA")
    assert str(rev.seq) == LEFT + sod.reverse_complement(SG) + RIGHT


def test_insert_does_not_mutate_input_vector():
    vec = _vector()
    before = str(vec.seq)
    sod.insert_sgrna(vec, STUFFER, SG, "sgA")
    assert str(vec.seq) == before


def test_features_inside_stuffer_removed_and_downstream_shifted():
    inside = SeqFeature(FeatureLocation(len(LEFT), len(LEFT) + len(STUFFER), 1),
                        type="misc_feature", qualifiers={"label": ["stuffer"]})
    down = SeqFeature(FeatureLocation(len(LEFT) + len(STUFFER) + 2, len(LEFT) + len(STUFFER) + 10, 1),
                      type="gene", qualifiers={"label": ["down"]})
    rec = sod.insert_sgrna(_vector([inside, down]), STUFFER, SG, "sgA")
    labels = [f.qualifiers["label"][0] for f in rec.features]
    assert "stuffer" not in labels and {"sgA", "down"} <= set(labels)
    d = next(f for f in rec.features if f.qualifiers["label"] == ["down"])
    assert str(d.extract(rec.seq)) == (LEFT + SG + RIGHT)[d.location.start:d.location.end]
    assert d.location.start == len(LEFT) + len(SG) + 2          # shifted by +3


@pytest.mark.xfail(strict=True, reason=(
    "Known bug: a feature spanning the insertion site keeps its old end coordinate "
    "when sgRNA and stuffer lengths differ. Remove this marker when fixed."))
def test_feature_spanning_stuffer_is_extended():
    span = SeqFeature(FeatureLocation(0, len(LEFT) + len(STUFFER) + 5, 1),
                      type="gene", qualifiers={"label": ["span"]})
    rec = sod.insert_sgrna(_vector([span]), STUFFER, SG, "sgA")
    f = next(x for x in rec.features if x.qualifiers["label"] == ["span"])
    assert f.location.end == len(LEFT) + len(SG) + 5


# --- I/O ------------------------------------------------------------------

def test_read_input_drops_blank_rows_and_bom(tmp_path):
    p = tmp_path / "in.csv"
    p.write_text("﻿sgRNA_id,sgRNA_seq\na," + SG + "\n,\n  ,  \n", encoding="utf-8")
    assert sod.read_input(p) == [{"sgRNA_id": "a", "sgRNA_seq": SG}]


def test_read_input_rejects_unknown_extension(tmp_path):
    p = tmp_path / "in.txt"
    p.write_text("x")
    with pytest.raises(ValueError, match="Unsupported"):
        sod.read_input(p)


def test_read_fasta_first_record_only(tmp_path):
    p = tmp_path / "s.fa"
    p.write_text(">one\nACGT\nAC\n>two\nGGGG\n")
    assert sod.read_fasta_sequence(p) == "ACGTAC"


def test_check_unique_ids():
    sod.check_unique_ids([{"sgRNA_id": "a"}, {"sgRNA_id": "b"}])
    with pytest.raises(ValueError, match="Duplicate sgRNA IDs"):
        sod.check_unique_ids([{"sgRNA_id": "a"}, {"sgRNA_id": "a"}])
    with pytest.raises(ValueError, match="not found"):
        sod.check_unique_ids([{"x": "a"}])

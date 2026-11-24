"""Unit tests for sgRNA cleaning/validation, IDT rows and file I/O."""

import pytest

import sgrna_oligo_designer as sod

SG = "GATTCTAAAACGGATTACCA"


def _site():
    """A minimal BsmBI-style site: CACC / GTTT overhangs."""
    return sod.CloningSite(enzyme="BsmBI", enzyme_site="CGTCTC", orientation="forward",
                           left_overhang="CACC", right_overhang="GTTT", replace_start=0,
                           replace_end=0, left_context="", right_context="", vector_length=1)


def test_reverse_complement():
    assert sod.reverse_complement("AACG") == "CGTT"


def test_default_oligos_match_the_classic_cacc_g_design():
    (row,) = sod.add_cloning_oligos([{"sgRNA_id": "a", "sgRNA_seq": SG}], _site())
    assert row["oligo_sense"] == "CACCG" + SG
    assert row["oligo_antisense"] == "AAAC" + sod.reverse_complement(SG) + "C"


def test_whitespace_and_case_cleaned():
    messy = " gattctaaaa\n cggattacca\t"
    (row,) = sod.add_cloning_oligos([{"sgRNA_id": "a", "sgRNA_seq": messy}], _site())
    assert row["sgRNA_seq"] == SG


@pytest.mark.parametrize("bad", ["", "   ", "GATTCTAAAACGGATTACCN", "GATTCTAAAACGGATTAC-A", "AUGCAUGCAUGCAUGCAUGC"])
def test_invalid_sgrna_raises(bad):
    with pytest.raises(ValueError):
        sod.add_cloning_oligos([{"sgRNA_id": "a", "sgRNA_seq": bad}], _site())


def test_missing_column_raises():
    with pytest.raises(ValueError, match="not found"):
        sod.add_cloning_oligos([{"id": "a", "seq": SG}], _site())


def test_custom_column_names():
    (row,) = sod.add_cloning_oligos([{"guide_id": "a", "guide": SG}], _site(),
                                    seq_col="guide", id_col="guide_id")
    assert row["oligo_sense"] == "CACCG" + SG


def test_idt_rows_interleaved():
    rows = sod.add_cloning_oligos([{"sgRNA_id": "a", "sgRNA_seq": SG},
                                   {"sgRNA_id": "b", "sgRNA_seq": SG}], _site())
    idt = sod.build_idt_rows(rows)
    assert [r["Name"] for r in idt] == ["a_s", "a_a", "b_s", "b_a"]
    assert {r["Scale"] for r in idt} == {"25nm"} and {r["Purification"] for r in idt} == {"STD"}


def test_check_unique_ids():
    sod.check_unique_ids([{"sgRNA_id": "a"}, {"sgRNA_id": "b"}])
    with pytest.raises(ValueError, match="Duplicate sgRNA IDs"):
        sod.check_unique_ids([{"sgRNA_id": "a"}, {"sgRNA_id": "a"}])
    with pytest.raises(ValueError, match="not found"):
        sod.check_unique_ids([{"x": "a"}])


def test_read_input_drops_blank_rows_and_bom(tmp_path):
    p = tmp_path / "in.csv"
    p.write_text("﻿sgRNA_id,sgRNA_seq\na," + SG + "\n,\n  ,  \n", encoding="utf-8")
    assert sod.read_input(p) == [{"sgRNA_id": "a", "sgRNA_seq": SG}]


def test_read_input_rejects_unknown_extension(tmp_path):
    p = tmp_path / "in.txt"
    p.write_text("x")
    with pytest.raises(ValueError, match="Unsupported"):
        sod.read_input(p)

"""Tests for the presentation helpers (tables, design template, download bundles)."""

import io
import zipfile
from pathlib import Path

import pytest

import report
import sgrna_oligo_designer as sod

ROOT = Path(__file__).resolve().parent.parent
SG = "GATTCTAAAACGGATTACCA"


@pytest.fixture(scope="module")
def ipusepr():
    vec = sod.read_vector_gb(ROOT / "PyShinyApp" / "vectors" / "ipUSEPR.gb")
    site = sod.find_cloning_site(vec, "BsmBI")
    rows = sod.add_cloning_oligos(
        [{"sgRNA_id": "a", "sgRNA_seq": SG}, {"sgRNA_id": "b", "sgRNA_seq": "AAAATTTTAAAACCCCGGGG"}], site)
    return vec, site, rows, sod.build_plasmid_records(rows, vec, site)


def test_gc_percent():
    assert report.gc_percent("GGCC") == 100.0 and report.gc_percent("ATAT") == 0.0
    assert report.gc_percent("") == 0.0 and report.gc_percent("GATC") == 50.0


def test_oligo_table_columns_and_notes(ipusepr):
    _, _, rows, _ = ipusepr
    t = report.oligo_table(rows)
    assert [r["sgRNA"] for r in t] == ["a", "b"]
    assert t[0]["Sense oligo (5'→3')"] == "CACCG" + SG and t[0]["Notes"] == "OK"
    assert t[0]["Length (nt)"] == "25 / 25"
    assert "TTTT" in t[1]["Notes"] and not t[1]["Notes"].startswith("sgRNA '")


def test_design_template_matches_the_classic_design(ipusepr):
    _, site, _, _ = ipusepr
    top, bottom = report.design_template(site, "always")
    assert top == "CACCG" + "N" * 20 and bottom == "AAAC" + "N" * 20 + "C"
    assert report.split_spacer(top) == ("CACCG", "N" * 20, "")
    assert report.design_template(site, "never")[0] == "CACC" + "N" * 20
    assert report.design_template(site, "if_needed") == report.design_template(site, "always")


def test_design_template_for_two_enzymes():
    rec = sod.read_vector_gb(ROOT / "tests" / "fixtures" / "pMCB306.gb")
    site = sod.find_cloning_site(rec, "BstXI", enzyme_b="BlpI")
    top, bottom = report.design_template(site, "never")
    assert report.split_spacer(top) == ("TTG", "N" * 20, "GTTTAAGAGC")
    assert bottom.startswith("TTAGCTCTTAAAC") and bottom.endswith("CAACAAG")


def test_oligos_csv_round_trip(ipusepr):
    import csv
    _, _, rows, _ = ipusepr
    got = list(csv.DictReader(io.StringIO(report.oligos_csv(rows).decode())))
    assert got[0]["sense_oligo"] == "CACCG" + SG and got[1]["sgRNA"] == "b"


def test_genbank_zip_has_one_named_file_per_sgrna(ipusepr):
    _, _, _, records = ipusepr
    z = zipfile.ZipFile(io.BytesIO(report.genbank_zip(records, "my vector")))
    assert sorted(z.namelist()) == ["my_vector_a.gb", "my_vector_b.gb"]
    assert z.read("my_vector_a.gb").decode().startswith("LOCUS       my_vector_a")


def test_bundle_zip_contents(ipusepr, tmp_path):
    _, _, rows, records = ipusepr
    z = zipfile.ZipFile(io.BytesIO(report.bundle_zip(
        sod.build_idt_rows(rows), rows, records, "ipUSEPR", "all plasmids", tmp_path / "o.xlsx")))
    assert set(z.namelist()) == {"IDT_orderForm.xlsx", "oligos.csv", "all_plasmids.gb",
                                 "plasmids/ipUSEPR_a.gb", "plasmids/ipUSEPR_b.gb"}
    assert z.read("IDT_orderForm.xlsx")[:2] == b"PK"           # a real xlsx (zip) file

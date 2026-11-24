"""
Golden tests: regenerate every file in example_outputs/ from the bundled vectors
and example sgRNA tables and require an exact match.

The golden plasmids were first produced by the old stuffer-based tool. The
enzyme-driven pipeline reproduces the same plasmid sequences; the only
intentional difference is the `source` feature, which now spans the whole
(longer or shorter) plasmid instead of keeping the empty vector's length.
"""

from pathlib import Path

import openpyxl
import pytest

import sgrna_oligo_designer as sod

ROOT    = Path(__file__).resolve().parent.parent
INPUTS  = ROOT / "example_inputs"
VECTORS = ROOT / "PyShinyApp" / "vectors"
GOLDEN  = ROOT / "example_outputs"

# (vector name, input table, vector .gb, enzyme, expected orientation)
CASES = [
    ("ipUSEPR",          "example_sgRNA_seqs.xlsx", "ipUSEPR.gb",          "BsmBI", "forward"),
    ("lentiCRISPRv2",    "example_sgRNA_seqs.csv",  "lentiCRISPRv2.gb",    "BsmBI", "forward"),
    ("LentiCRISPRv2GFP", "example_sgRNA_seqs.tsv",  "LentiCRISPRv2GFP.gb", "BsmBI", "reverse"),
]
IDS = [c[0] for c in CASES]


def _run(table, vector, enzyme):
    vec  = sod.read_vector_gb(VECTORS / vector)
    site = sod.find_cloning_site(vec, enzyme)
    rows = sod.add_cloning_oligos(sod.read_input(INPUTS / "sgRNA_seqs" / table), site)
    return rows, vec, site, sod.build_plasmid_records(rows, vec, site)


def _sheet(path):
    wb = openpyxl.load_workbook(path, data_only=True)
    data = [tuple(r) for r in wb.active.iter_rows(values_only=True)]
    wb.close()
    return data


@pytest.mark.parametrize("name,table,vector,enzyme,orientation", CASES, ids=IDS)
def test_per_sgrna_genbank_files_match_golden(name, table, vector, enzyme, orientation, tmp_path):
    _, _, site, records = _run(table, vector, enzyme)
    assert site.orientation == orientation
    sod.write_gb_files(records, tmp_path, vector_name=name)

    written  = sorted(p.name for p in tmp_path.glob("*.gb"))
    expected = sorted(p.name for p in (GOLDEN / "plasmids").glob(f"{name}_sg_*.gb"))
    assert written == expected and len(written) == 4
    for fname in written:
        assert (tmp_path / fname).read_text() == (GOLDEN / "plasmids" / fname).read_text(), fname


def test_combined_genbank_matches_golden(tmp_path):
    name, table, vector, enzyme, _ = CASES[0]
    _, _, _, records = _run(table, vector, enzyme)
    sod.write_gb_files(records, tmp_path, vector_name=name, multi_gb=True,
                       combined_stem="ipUSEPR_plasmids_combined")
    got = (tmp_path / "ipUSEPR_plasmids_combined.gb").read_text()
    assert got == (GOLDEN / "plasmids" / "ipUSEPR_plasmids_combined.gb").read_text()


def test_idt_order_form_matches_golden(tmp_path):
    name, table, vector, enzyme, _ = CASES[0]
    rows, _, _, _ = _run(table, vector, enzyme)
    out = tmp_path / "order.xlsx"
    sod.write_idt_excel(sod.build_idt_rows(rows), out)
    assert _sheet(out) == _sheet(GOLDEN / "example_sgRNA_seqs_IDT_orderForm.xlsx")


@pytest.mark.parametrize("name,table,vector,enzyme,orientation", CASES, ids=IDS)
def test_all_example_vectors_give_the_same_overhangs_and_no_warnings(
        name, table, vector, enzyme, orientation):
    rows, _, site, _ = _run(table, vector, enzyme)
    assert (site.left_overhang, site.right_overhang) == ("CACC", "GTTT")
    assert not [w for r in rows for w in r["warnings"]]
    for r in rows:
        assert r["oligo_sense"] == "CACCG" + r["sgRNA_seq"]
        assert r["oligo_antisense"] == "AAAC" + sod.reverse_complement(r["sgRNA_seq"]) + "C"


@pytest.mark.parametrize("table", ["example_sgRNA_seqs.xlsx", "example_sgRNA_seqs.csv",
                                   "example_sgRNA_seqs.tsv"])
def test_all_input_formats_give_same_rows(table):
    ref = sod.read_input(INPUTS / "sgRNA_seqs" / "example_sgRNA_seqs.xlsx")
    got = sod.read_input(INPUTS / "sgRNA_seqs" / table)
    key = lambda rows: [(r["sgRNA_id"], r["sgRNA_seq"]) for r in rows]
    assert key(got) == key(ref)


@pytest.mark.parametrize("name,table,vector,enzyme,orientation", CASES, ids=IDS)
def test_plasmid_structure(name, table, vector, enzyme, orientation):
    """Independent of golden files: spacer present, length right, no stale features."""
    rows, vec, site, records = _run(table, vector, enzyme)
    for row, rec in zip(rows, records):
        sg = row["sgRNA_seq"]
        assert len(rec.seq) == len(vec.seq) - site.removed_length + len(row["insert_core"])
        feats = [f for f in rec.features if f.qualifiers.get("label") == [row["sgRNA_id"]]]
        assert len(feats) == 1
        assert str(feats[0].extract(rec.seq)) == sg          # strand-aware read-back
        # every 'source' feature spans the whole plasmid
        for f in rec.features:
            if f.type == "source" and len(f.location.parts) == 1:
                assert (int(f.location.start), int(f.location.end)) == (0, len(rec.seq))


def test_bundled_example_table_matches_repo_example():
    """The deploy dir is self-contained, so it carries its own copy of the example table."""
    bundled = ROOT / "PyShinyApp" / "example_inputs" / "example_sgRNA_seqs.xlsx"
    assert _sheet(bundled) == _sheet(INPUTS / "sgRNA_seqs" / "example_sgRNA_seqs.xlsx")


def test_app_imports_and_builds():
    """Smoke test: the Shiny app module imports and exposes an App."""
    import app
    from shiny import App
    assert isinstance(app.app, App)

"""
Golden tests: regenerate every file in example_outputs/ from example_inputs/
and require an exact match. These guard the uv migration, restructure,
dependency upgrades and the enzyme-based redesign (for the existing vectors).
"""

from pathlib import Path

import openpyxl
import pytest

import sgrna_oligo_designer as sod

ROOT = Path(__file__).resolve().parent.parent
INPUTS = ROOT / "example_inputs"
GOLDEN = ROOT / "example_outputs"

# (vector name, input table, vector .gb, stuffer: raw string or FASTA path)
CASES = [
    ("ipUSEPR", "example_sgRNA_seqs.xlsx", "ipUSEPR.gb", "GAGACGTGACCGTCTCT"),
    ("lentiCRISPRv2", "example_sgRNA_seqs.csv", "lentiCRISPRv2.gb",
     INPUTS / "stuffer_seqs" / "lentiCRISPRv2_stuffer.fa"),
    ("LentiCRISPRv2GFP", "example_sgRNA_seqs.tsv", "LentiCRISPRv2GFP.gb",
     INPUTS / "stuffer_seqs" / "LentiCRISPRv2GFP_stuffer.fa"),
]
IDS = [c[0] for c in CASES]


def _stuffer(value) -> str:
    return sod.read_fasta_sequence(Path(value)) if isinstance(value, Path) else value


def _run(table, vector, stuffer):
    rows = sod.read_input(INPUTS / "sgRNA_seqs" / table)
    rows = sod.add_cloning_oligos(rows)
    vec = sod.read_vector_gb(INPUTS / "vector_seqs" / vector)
    records = sod.build_plasmid_records(rows, vec, _stuffer(stuffer))
    return rows, vec, records


def _sheet(path):
    wb = openpyxl.load_workbook(path, data_only=True)
    ws = wb.active
    data = [tuple(r) for r in ws.iter_rows(values_only=True)]
    wb.close()
    return data


@pytest.mark.parametrize("name,table,vector,stuffer", CASES, ids=IDS)
def test_per_sgrna_genbank_files_match_golden(name, table, vector, stuffer, tmp_path):
    _, _, records = _run(table, vector, stuffer)
    sod.write_gb_files(records, tmp_path, vector_name=name)

    written = sorted(p.name for p in tmp_path.glob("*.gb"))
    expected = sorted(p.name for p in (GOLDEN / "plasmids").glob(f"{name}_sg_*.gb"))
    assert written == expected and len(written) == 4
    for fname in written:
        assert (tmp_path / fname).read_text() == (GOLDEN / "plasmids" / fname).read_text(), fname


def test_combined_genbank_matches_golden(tmp_path):
    name, table, vector, stuffer = CASES[0]
    _, _, records = _run(table, vector, stuffer)
    sod.write_gb_files(records, tmp_path, vector_name=name, multi_gb=True,
                       combined_stem="ipUSEPR_plasmids_combined")
    got = (tmp_path / "ipUSEPR_plasmids_combined.gb").read_text()
    assert got == (GOLDEN / "plasmids" / "ipUSEPR_plasmids_combined.gb").read_text()


def test_idt_order_form_matches_golden(tmp_path):
    rows = sod.add_cloning_oligos(
        sod.read_input(INPUTS / "sgRNA_seqs" / "example_sgRNA_seqs.xlsx"))
    out = tmp_path / "order.xlsx"
    sod.write_idt_excel(sod.build_idt_rows(rows), out)
    assert _sheet(out) == _sheet(GOLDEN / "example_sgRNA_seqs_IDT_orderForm.xlsx")


@pytest.mark.parametrize("table", ["example_sgRNA_seqs.xlsx", "example_sgRNA_seqs.csv",
                                   "example_sgRNA_seqs.tsv"])
def test_all_input_formats_give_same_rows(table):
    ref = sod.read_input(INPUTS / "sgRNA_seqs" / "example_sgRNA_seqs.xlsx")
    got = sod.read_input(INPUTS / "sgRNA_seqs" / table)
    key = lambda rows: [(r["sgRNA_id"], r["sgRNA_seq"]) for r in rows]
    assert key(got) == key(ref)


@pytest.mark.parametrize("name,table,vector,stuffer", CASES, ids=IDS)
def test_plasmid_structure(name, table, vector, stuffer):
    """Independent of golden files: sgRNA present, length and strand consistent."""
    rows, vec, records = _run(table, vector, stuffer)
    stuffer_len = len(_stuffer(stuffer))
    for row, rec in zip(rows, records):
        sg = row["sgRNA_seq"]
        assert len(rec.seq) == len(vec.seq) - stuffer_len + len(sg)
        feat = [f for f in rec.features if f.qualifiers.get("label") == [row["sgRNA_id"]]]
        assert len(feat) == 1
        f = feat[0]
        got = str(f.extract(rec.seq))   # strand-aware: always reads as the sgRNA
        # Feature on -1 strand extracts as reverse complement of stored bases == sgRNA.
        assert got == sg
        assert [x for x in rec.features if x.location.start == f.location.start
                and x.location.end == f.location.end] == [f]


def test_root_and_shiny_copies_are_identical():
    """Temporary guard until the restructure leaves a single copy."""
    a = (ROOT / "sgrna_oligo_designer.py").read_bytes()
    b = (ROOT / "PyShinyApp" / "sgrna_oligo_designer.py").read_bytes()
    assert a == b

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


def test_validate_ids_accepts_normal_ids_and_strips_surrounding_whitespace():
    rows = [{"sgRNA_id": " sg_Luc-474.v2 "}, {"sgRNA_id": "b"}]
    sod.validate_ids(rows)
    assert rows[0]["sgRNA_id"] == "sg_Luc-474.v2"


@pytest.mark.parametrize("bad", ["my guide", "=HYPERLINK(1)", "+cmd", "-1", ".hidden", "a/b", "a;b", "", "  ", "sg\u00e9"])
def test_validate_ids_rejects_unsafe_ids(bad):
    with pytest.raises(ValueError, match="Invalid sgRNA ID"):
        sod.validate_ids([{"sgRNA_id": "ok"}, {"sgRNA_id": bad}])


def test_validate_ids_duplicates_and_missing_column():
    with pytest.raises(ValueError, match="Duplicate sgRNA IDs"):
        sod.validate_ids([{"sgRNA_id": "a"}, {"sgRNA_id": "a"}])
    with pytest.raises(ValueError, match="ignoring upper/lower case"):         # SG1 / sg1 collide as files
        sod.validate_ids([{"sgRNA_id": "SG1"}, {"sgRNA_id": "sg1"}])
    with pytest.raises(ValueError, match="not found"):
        sod.validate_ids([{"x": "a"}])


@pytest.mark.parametrize("raw,expected", [
    ("my vector", "my_vector"), ("../../evil", "evil"), ("/abs/path/x", "abs_path_x"),
    ("..", "unnamed"), ("", "unnamed"), ("a b/c\\d:e", "a_b_c_d_e"), ("ok-name_1.2", "ok-name_1.2"),
])
def test_sanitize_name(raw, expected):
    assert sod.sanitize_name(raw) == expected
    assert "/" not in sod.sanitize_name(raw) and "\\" not in sod.sanitize_name(raw)


def test_idt_excel_never_writes_a_formula(tmp_path):
    import openpyxl
    rows = [{"Name": "=HYPERLINK(\"http://x\",\"a\")", "Sequence": "ACGT", "Scale": "25nm", "Purification": "STD"}]
    out = tmp_path / "o.xlsx"
    sod.write_idt_excel(rows, out)
    cell = openpyxl.load_workbook(out).active["A2"]
    assert cell.data_type == "s" and cell.value.startswith("=")


def _plasmids():
    from pathlib import Path
    vec = sod.read_vector_gb(Path(__file__).resolve().parent.parent / "PyShinyApp" / "vectors" / "ipUSEPR.gb")
    site = sod.find_cloning_site(vec, "BsmBI")
    rows = sod.add_cloning_oligos([{"sgRNA_id": "a", "sgRNA_seq": SG}], site)
    return sod.build_plasmid_records(rows, vec, site)


@pytest.mark.parametrize("stem", ["../../escaped", "/tmp/escaped_abs", "a/b/c", "..", ""])
def test_combined_filename_cannot_leave_the_output_directory(tmp_path, stem):
    out = tmp_path / "out"
    sod.write_gb_files(_plasmids(), out, "v", multi_gb=True, combined_stem=stem)
    written = list(tmp_path.rglob("*.gb"))
    assert len(written) == 1 and written[0].parent == out


def test_vector_name_with_spaces_no_longer_breaks_genbank_output(tmp_path):
    sod.write_gb_files(_plasmids(), tmp_path, "my vector v2", multi_gb=True)
    text = next(tmp_path.glob("*.gb")).read_text()
    assert "my_vector_v2_a" in text.splitlines()[0]


def test_read_input_drops_blank_rows_and_bom(tmp_path):
    p = tmp_path / "in.csv"
    p.write_text("﻿sgRNA_id,sgRNA_seq\na," + SG + "\n,\n  ,  \n", encoding="utf-8")
    assert sod.read_input(p) == [{"sgRNA_id": "a", "sgRNA_seq": SG}]


def test_read_input_rejects_unknown_extension(tmp_path):
    p = tmp_path / "in.txt"
    p.write_text("x")
    with pytest.raises(ValueError, match="Unsupported"):
        sod.read_input(p)


# --- read_input edge cases -----------------------------------------------------------

def _xlsx(path, rows):
    import openpyxl
    wb = openpyxl.Workbook()
    for r in rows:
        wb.active.append(r)
    wb.save(path)


def test_empty_spreadsheet_gives_a_clear_error(tmp_path):
    import openpyxl
    p = tmp_path / "empty.xlsx"
    openpyxl.Workbook().save(p)
    with pytest.raises(ValueError, match="No sgRNA rows|empty"):
        sod.read_input(p)


def test_header_only_and_blank_csv_give_a_clear_error(tmp_path):
    p = tmp_path / "h.csv"
    p.write_text("sgRNA_id,sgRNA_seq\n,\n")
    with pytest.raises(ValueError, match="No sgRNA rows"):
        sod.read_input(p)


def test_blank_header_columns_are_dropped(tmp_path):
    p = tmp_path / "t.xlsx"
    _xlsx(p, [["sgRNA_id", "sgRNA_seq", None], ["a", SG, "stray"]])
    assert sod.read_input(p) == [{"sgRNA_id": "a", "sgRNA_seq": SG}]


def test_row_cap(tmp_path):
    p = tmp_path / "big.csv"
    p.write_text("sgRNA_id,sgRNA_seq\n" + "".join(f"s{i},{SG}\n" for i in range(sod.MAX_SGRNAS + 1)))
    with pytest.raises(ValueError, match="up to 1000"):
        sod.read_input(p)


def test_xls_is_not_claimed_to_be_supported(tmp_path):
    p = tmp_path / "old.xls"
    p.write_bytes(b"x")
    with pytest.raises(ValueError, match="Unsupported"):
        sod.read_input(p)

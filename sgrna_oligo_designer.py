"""
sgrna_oligo_designer.py
------------------------
Generates annealing oligo pairs for sgRNA cloning from a spreadsheet,
and produces a ready-to-upload IDT bulk order Excel file by default.

Dependencies: openpyxl only (stdlib csv/pathlib handle everything else).

Default parameters target the ipUSEPR vector (BsmBI / U6 promoter):
  - Top oligo  : 5'-CACCG-{sgRNA_seq}-3'
  - Bottom oligo: 5'-AAAC-{rc(sgRNA_seq)}C-3'

Usage (CLI):
    # Default: generates both CSV and IDT Excel
    python sgrna_oligo_designer.py input.xlsx

    # Custom output directory (both files land there)
    python sgrna_oligo_designer.py input.xlsx --output /path/to/results/

    # CSV only, skip IDT file
    python sgrna_oligo_designer.py input.xlsx --no-idt

    # Custom vector overrides
    python sgrna_oligo_designer.py input.xlsx --top-prefix CACC --bot-prefix AAAC --bot-suffix ""

Usage (Python):
    from sgrna_oligo_designer import add_cloning_oligos, write_idt_excel

    rows = add_cloning_oligos(rows)                     # ipUSEPR defaults
    write_idt_excel(rows, Path("my_order_IDT.xlsx"))    # IDT-ready Excel
"""

import argparse
import csv
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# IDT template constants (matches IDT bulk upload spec exactly)
# ---------------------------------------------------------------------------

IDT_COLUMNS    = ["Name", "Sequence", "Scale", "Purification"]
IDT_SCALE      = "25nm"
IDT_PURIF      = "STD"
IDT_COL_WIDTHS = {"A": 12.54, "B": 41.54, "C": 10.18, "D": 18.82}


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence (ACGT only)."""
    return seq.translate(_COMPLEMENT)[::-1]


def add_cloning_oligos(
    rows: list[dict],
    seq_col: str = "sgRNA_seq",
    top_prefix: str = "CACCG",
    bot_prefix: str = "AAAC",
    bot_suffix: str = "C",
) -> list[dict]:
    """
    Add 'oligo_top' and 'oligo_bottom' keys to each row dict.

    Parameters
    ----------
    rows : list[dict]
        Input rows, each containing at minimum the sgRNA sequence column.
    seq_col : str
        Key holding the raw sgRNA sequences (20-mer, no PAM).
    top_prefix : str
        5' prefix appended to the sense strand.
        ipUSEPR default: 'CACCG'  (CACC overhang + G transcriptional start)
    bot_prefix : str
        5' prefix appended to the antisense strand.
        ipUSEPR default: 'AAAC'   (complement of GTTT BsmBI overhang)
    bot_suffix : str
        3' suffix appended to the antisense strand.
        ipUSEPR default: 'C'      (complement of the extra G on top strand)

    Returns
    -------
    list[dict]
        Same rows with 'oligo_top' and 'oligo_bottom' added in-place.
    """
    if rows and seq_col not in rows[0]:
        raise ValueError(
            f"Column '{seq_col}' not found. Available: {list(rows[0].keys())}"
        )

    for row in rows:
        seq = row[seq_col]
        row["oligo_top"]    = top_prefix + seq
        row["oligo_bottom"] = bot_prefix + reverse_complement(seq) + bot_suffix

    return rows


def build_idt_rows(
    rows: list[dict],
    id_col: str = "sgRNA_id",
) -> list[dict]:
    """
    Reshape oligo rows into IDT bulk upload format.

    Each sgRNA produces two interleaved rows (top then bottom):
        <sgRNA_id>_top    | oligo_top    | 25nm | STD
        <sgRNA_id>_bottom | oligo_bottom | 25nm | STD

    Parameters
    ----------
    rows : list[dict]
        Rows with sgRNA_id, oligo_top, and oligo_bottom keys.
    id_col : str
        Key holding the sgRNA identifiers.

    Returns
    -------
    list[dict]
        IDT-formatted rows with keys: Name, Sequence, Scale, Purification.
    """
    idt_rows = []
    for r in rows:
        idt_rows.append({
            "Name":         f"{r[id_col]}_top",
            "Sequence":     r["oligo_top"],
            "Scale":        IDT_SCALE,
            "Purification": IDT_PURIF,
        })
        idt_rows.append({
            "Name":         f"{r[id_col]}_bottom",
            "Sequence":     r["oligo_bottom"],
            "Scale":        IDT_SCALE,
            "Purification": IDT_PURIF,
        })
    return idt_rows


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def read_input(path: Path) -> list[dict]:
    """
    Read an input spreadsheet or delimited file into a list of row dicts.
    Supports .xlsx / .xlsm / .xls (via openpyxl), .csv, and .tsv.
    """
    suffix = path.suffix.lower()

    if suffix in {".xlsx", ".xlsm", ".xls"}:
        import openpyxl
        wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
        ws = wb.active
        rows_iter = ws.iter_rows(values_only=True)
        headers = [str(h) for h in next(rows_iter)]
        rows = [dict(zip(headers, (str(v) if v is not None else "" for v in row)))
                for row in rows_iter]
        wb.close()
        return rows

    elif suffix == ".csv":
        with open(path, newline="", encoding="utf-8-sig") as fh:
            return list(csv.DictReader(fh))

    elif suffix == ".tsv":
        with open(path, newline="", encoding="utf-8-sig") as fh:
            return list(csv.DictReader(fh, delimiter="\t"))

    else:
        raise ValueError(f"Unsupported file type: '{suffix}'. Use .xlsx, .csv, or .tsv")


def write_csv_output(rows: list[dict], path: Path) -> None:
    """Write full oligo rows to a flat CSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        print(f"Saved → {path} (empty)")
        return
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved → {path}")


def write_idt_excel(idt_rows: list[dict], path: Path) -> None:
    """
    Write IDT-formatted rows to Excel, replicating the official IDT bulk
    upload template (bold Calibri headers, IDT column widths, no index).
    """
    import openpyxl
    from openpyxl.styles import Font

    path.parent.mkdir(parents=True, exist_ok=True)

    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Sheet1"

    # Header row — bold Calibri, matching IDT template exactly
    header_font = Font(name="Calibri", bold=True)
    for col_idx, col_name in enumerate(IDT_COLUMNS, start=1):
        cell = ws.cell(row=1, column=col_idx, value=col_name)
        cell.font = header_font

    # Data rows
    data_font = Font(name="Calibri")
    for row_idx, row in enumerate(idt_rows, start=2):
        for col_idx, col_name in enumerate(IDT_COLUMNS, start=1):
            cell = ws.cell(row=row_idx, column=col_idx, value=row[col_name])
            cell.font = data_font

    # Column widths — preserved from original IDT template
    for col_letter, width in IDT_COL_WIDTHS.items():
        ws.column_dimensions[col_letter].width = width

    wb.save(path)
    print(f"Saved → {path}")


def resolve_output_path(args_output: Path | None, input_path: Path, suffix: str) -> Path:
    """
    Resolve an output path given the --output argument, input path, and
    a filename suffix (e.g. '_oligos.csv' or '_oligos_IDT.xlsx').
    Handles three cases:
      - No --output:           same directory as input, default filename
      - Directory / trailing slash: that directory, default filename
      - Explicit file path:    used as-is (suffix ignored for CSV, swapped for IDT)
    """
    default_name = input_path.stem + suffix
    if args_output is None:
        return input_path.with_name(default_name)
    elif args_output.is_dir() or str(args_output).endswith(("/", "\\")):
        return args_output / default_name
    elif args_output.suffix.lower() not in {".csv", ".xlsx"}:
        return args_output / default_name
    else:
        if suffix.endswith(".xlsx"):
            return args_output.with_name(args_output.stem + "_IDT.xlsx")
        return args_output


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="sgrna_oligo_designer",
        description=(
            "Generate BsmBI/Golden Gate annealing oligos from an sgRNA spreadsheet. "
            "Produces a flat CSV and an IDT bulk-order Excel file by default."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "input", type=Path,
        help="Input file (.xlsx, .csv, or .tsv)",
    )
    p.add_argument(
        "--output", "-o", type=Path, default=None,
        help=(
            "Output directory or full CSV path. "
            "IDT Excel is written alongside with '_IDT.xlsx' suffix. "
            "Default: same directory as input."
        ),
    )
    p.add_argument(
        "--no-idt", action="store_true",
        help="Skip IDT Excel output and generate CSV only.",
    )
    p.add_argument("--seq-col",    default="sgRNA_seq", help="Column name for sgRNA sequences")
    p.add_argument("--id-col",     default="sgRNA_id",  help="Column name for sgRNA identifiers")
    p.add_argument("--top-prefix", default="CACCG",     help="5' prefix for top (sense) oligo")
    p.add_argument("--bot-prefix", default="AAAC",      help="5' prefix for bottom (antisense) oligo")
    p.add_argument("--bot-suffix", default="C",         help="3' suffix for bottom (antisense) oligo")
    return p


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    input_path: Path = args.input
    if not input_path.exists():
        print(f"Error: file not found — {input_path}", file=sys.stderr)
        sys.exit(1)

    csv_path = resolve_output_path(args.output, input_path, "_oligos.csv")
    idt_path = resolve_output_path(args.output, input_path, "_oligos_IDT.xlsx")

    # Read, process
    rows = read_input(input_path)
    rows = add_cloning_oligos(
        rows,
        seq_col=args.seq_col,
        top_prefix=args.top_prefix,
        bot_prefix=args.bot_prefix,
        bot_suffix=args.bot_suffix,
    )

    # Write CSV (always)
    write_csv_output(rows, csv_path)

    # Write IDT Excel (default, unless --no-idt)
    if not args.no_idt:
        idt_rows = build_idt_rows(rows, id_col=args.id_col)
        write_idt_excel(idt_rows, idt_path)

    # Summary
    n = len(rows)
    print(f"Processed {n} sgRNA sequences → {n * 2} oligos to order ({n} top + {n} bottom)")


if __name__ == "__main__":
    main()

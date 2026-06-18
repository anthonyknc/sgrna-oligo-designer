"""
sgrna_oligo_designer.py
------------------------
Generates annealing oligo pairs for sgRNA cloning from a spreadsheet,
writes a ready-to-upload IDT bulk order Excel file, and optionally
produces annotated full-plasmid GenBank files for sequencing verification.

Dependencies: openpyxl (always), biopython (only for .gb output).

Default parameters target the ipUSEPR vector (BsmBI / U6 promoter):
  - Top oligo  : 5'-CACCG-{sgRNA_seq}-3'
  - Bottom oligo: 5'-AAAC-{rc(sgRNA_seq)}C-3'

Usage:
    # IDT Excel only (no biopython required)
    python sgrna_oligo_designer.py input.xlsx

    # IDT Excel + one .gb per sgRNA (default)
    python sgrna_oligo_designer.py input.xlsx --vector empty_vector.gb

    # IDT Excel + single combined .gb
    python sgrna_oligo_designer.py input.xlsx --vector empty_vector.gb --multi-gb

    # Custom stuffer sequence in the vector
    python sgrna_oligo_designer.py input.xlsx --vector empty_vector.gb --stuffer ATCGATCGATCGATCGATCG

    # Custom output locations
    python sgrna_oligo_designer.py input.xlsx --output results/ --vector empty_vector.gb --gb-output results/plasmids/

    # Custom vector overrides
    python sgrna_oligo_designer.py input.xlsx --top-prefix CACC --bot-prefix AAAC --bot-suffix ""
"""

import argparse
import copy
import csv
import sys
from collections import Counter
from pathlib import Path


# ---------------------------------------------------------------------------
# IDT template constants (matches IDT bulk upload spec exactly)
# ---------------------------------------------------------------------------

IDT_COLUMNS    = ["Name", "Sequence", "Scale", "Purification"]
IDT_SCALE      = "25nm"
IDT_PURIF      = "STD"
IDT_COL_WIDTHS = {"A": 12.54, "B": 41.54, "C": 10.18, "D": 18.82}

DEFAULT_STUFFER = "NNNNNNNNNNNNNNNNNNNN"   # 20 N's — common vector placeholder


# ---------------------------------------------------------------------------
# Core logic — oligo design
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
            "Name":         f"{r[id_col]}_s",
            "Sequence":     r["oligo_top"],
            "Scale":        IDT_SCALE,
            "Purification": IDT_PURIF,
        })
        idt_rows.append({
            "Name":         f"{r[id_col]}_a",
            "Sequence":     r["oligo_bottom"],
            "Scale":        IDT_SCALE,
            "Purification": IDT_PURIF,
        })
    return idt_rows


# ---------------------------------------------------------------------------
# Core logic — GenBank plasmid construction
# ---------------------------------------------------------------------------

def read_vector_gb(path: Path):
    """
    Parse an empty vector GenBank file.

    Parameters
    ----------
    path : Path
        Path to the .gb / .gbk file.

    Returns
    -------
    Bio.SeqRecord.SeqRecord
        Parsed vector record with sequence and features.
    """
    from Bio import SeqIO
    return SeqIO.read(path, "genbank")


def _find_all(seq: str, pattern: str) -> list[int]:
    """Return all 0-based start positions of pattern in seq."""
    positions, start = [], 0
    while (pos := seq.find(pattern, start)) != -1:
        positions.append(pos)
        start = pos + 1
    return positions


def _validate_stuffer_sequence(stuffer_seq: str) -> None:
    """
    Validate that the stuffer is non-empty and contains only ACGT bases.

    Raises
    ------
    ValueError
        If the stuffer is empty or contains non-ACGT characters.
    """
    if not stuffer_seq:
        raise ValueError(
            "Stuffer sequence is empty. "
            "Check --stuffer provides a valid sequence or non-empty FASTA file."
        )
    invalid = set(stuffer_seq.upper()) - set("ACGT")
    if invalid:
        raise ValueError(
            f"Stuffer contains non-ACGT characters: {sorted(invalid)}. "
            "Only standard ACGT bases are supported. "
            "Check your --stuffer input or FASTA file."
        )


def validate_stuffer(vector_record, stuffer_seq: str) -> tuple[int, int]:
    """
    Locate the stuffer sequence in the vector across both strands.
    Must appear exactly once across forward and reverse strands combined.

    Parameters
    ----------
    vector_record : Bio.SeqRecord.SeqRecord
        Parsed empty vector.
    stuffer_seq : str
        Placeholder sequence to find and replace (ACGT only).

    Returns
    -------
    tuple[int, int]
        (start_position, strand) where start_position is the 0-based position
        on the forward strand, and strand is +1 (forward) or -1 (reverse).

    Raises
    ------
    ValueError
        If the stuffer is empty, non-ACGT, palindromic, absent, or non-unique.
    """
    _validate_stuffer_sequence(stuffer_seq)

    fwd_seq     = str(vector_record.seq).upper()
    stuffer_fwd = stuffer_seq.upper()
    stuffer_rev = reverse_complement(stuffer_fwd)

    # Palindrome check — rc equals itself, strand is indeterminate
    if stuffer_fwd == stuffer_rev:
        raise ValueError(
            f"Stuffer sequence is palindromic (its reverse complement equals itself). "
            "Strand cannot be determined. "
            "Please provide a longer, asymmetric stuffer sequence."
        )

    fwd_hits = _find_all(fwd_seq, stuffer_fwd)
    rev_hits = _find_all(fwd_seq, stuffer_rev)
    total    = len(fwd_hits) + len(rev_hits)

    if total == 0:
        raise ValueError(
            f"Stuffer not found on either strand of vector '{vector_record.name}'. "
            "Check --stuffer matches the placeholder in your .gb file. "
            "Note: stuffers that span the circular origin are not supported."
        )
    if total > 1:
        raise ValueError(
            f"Stuffer found {total} times across both strands "
            f"(forward: {len(fwd_hits)}, reverse: {len(rev_hits)}). "
            "Stuffer must be unique across both strands."
        )

    if fwd_hits:
        return fwd_hits[0], +1
    else:
        return rev_hits[0], -1


def _shift_feature(feature, delta: int):
    """
    Return a copy of a BioPython SeqFeature with coordinates shifted by delta.
    Handles both simple FeatureLocation and CompoundLocation (e.g. join()).
    """
    from Bio.SeqFeature import FeatureLocation, CompoundLocation

    new_feat = copy.deepcopy(feature)
    loc      = feature.location

    if isinstance(loc, CompoundLocation):
        new_feat.location = CompoundLocation(
            [FeatureLocation(p.start + delta, p.end + delta, p.strand)
             for p in loc.parts],
            operator=loc.operator,
        )
    else:
        new_feat.location = FeatureLocation(
            loc.start + delta, loc.end + delta, loc.strand
        )
    return new_feat


def insert_sgrna(
    vector_record,
    stuffer_seq: str,
    sgrna_seq: str,
    sgrna_id: str,
):
    """
    Replace the stuffer with the raw 20-mer sgRNA sequence and annotate it.

    The stuffer is replaced in-place; all vector features downstream of the
    insertion site are coordinate-shifted if stuffer and sgRNA lengths differ
    (typically they are identical at 20 bp, so no shift is needed).
    A new misc_RNA feature labelled with sgrna_id is added at the insertion site.

    Parameters
    ----------
    vector_record : Bio.SeqRecord.SeqRecord
        Parsed empty vector (not modified in-place).
    stuffer_seq : str
        Placeholder sequence to replace.
    sgrna_seq : str
        Raw 20-mer sgRNA sequence (no PAM, no overhangs).
    sgrna_id : str
        Identifier for the sgRNA; used as the feature label and record ID.

    Returns
    -------
    Bio.SeqRecord.SeqRecord
        New record with sgRNA inserted and annotated.
    """
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, FeatureLocation

    insert_start, strand = validate_stuffer(vector_record, stuffer_seq)
    insert_end            = insert_start + len(stuffer_seq)

    # Insert the sgRNA in the correct orientation relative to the stuffer strand
    insert_seq = sgrna_seq.upper() if strand == +1 else reverse_complement(sgrna_seq.upper())
    delta      = len(insert_seq) - len(stuffer_seq)

    # Build new sequence
    old_seq = str(vector_record.seq)
    new_seq = old_seq[:insert_start] + insert_seq + old_seq[insert_end:]

    # Deep copy record, update sequence and metadata
    new_record             = copy.deepcopy(vector_record)
    new_record.seq         = Seq(new_seq)
    new_record.id          = sgrna_id
    new_record.name        = sgrna_id[:16]   # LOCUS name field: max 16 chars
    new_record.description = (
        f"{sgrna_id} cloned into {vector_record.name.strip()}"
    )

    # Shift downstream features if insertion changes total length
    if delta != 0:
        new_record.features = [
            _shift_feature(f, delta) if f.location.start >= insert_end else f
            for f in new_record.features
        ]

    # Annotate the inserted sgRNA on the correct strand
    new_record.features.append(
        SeqFeature(
            FeatureLocation(insert_start, insert_start + len(insert_seq), strand=strand),
            type="misc_feature",
            qualifiers={
                "label": [sgrna_id],
                "note":  [f"sgRNA sequence: {sgrna_seq}"],
                "gene":  [sgrna_id],
            },
        )
    )

    return new_record


def build_plasmid_records(
    rows: list[dict],
    vector_record,
    stuffer_seq: str,
    seq_col: str = "sgRNA_seq",
    id_col: str  = "sgRNA_id",
) -> list:
    """
    Generate one annotated SeqRecord per sgRNA row.

    Parameters
    ----------
    rows : list[dict]
        Rows with sgRNA_id and sgRNA_seq keys.
    vector_record : Bio.SeqRecord.SeqRecord
        Parsed empty vector.
    stuffer_seq : str
        Placeholder sequence to replace in the vector.
    seq_col : str
        Key holding raw sgRNA sequences.
    id_col : str
        Key holding sgRNA identifiers.

    Returns
    -------
    list[Bio.SeqRecord.SeqRecord]
        One record per sgRNA with the stuffer replaced and feature annotated.
    """
    return [
        insert_sgrna(vector_record, stuffer_seq, row[seq_col], row[id_col])
        for row in rows
    ]


def write_gb_files(
    records: list,
    output_dir: Path,
    vector_name: str,
    multi_gb: bool = False,
    combined_stem: str = "combined_plasmids",
) -> None:
    """
    Write GenBank files — one per sgRNA (default) or a single combined file.

    Parameters
    ----------
    records : list[Bio.SeqRecord.SeqRecord]
        Annotated plasmid records to write.
    output_dir : Path
        Directory for output .gb files.
    vector_name : str
        Vector name appended to each output filename.
        Produces: {sgRNA_id}_{vector_name}.gb
    multi_gb : bool
        If True, write all records into a single combined .gb file.
    combined_stem : str
        Filename stem for the combined file (used when multi_gb=True).
    """
    from Bio import SeqIO

    output_dir.mkdir(parents=True, exist_ok=True)

    if multi_gb:
        out_path = output_dir / f"{combined_stem}.gb"
        with open(out_path, "w") as fh:
            SeqIO.write(records, fh, "genbank")
        print(f"Saved → {out_path}")
    else:
        for record in records:
            out_path = output_dir / f"{vector_name}_{record.id}.gb"
            with open(out_path, "w") as fh:
                SeqIO.write([record], fh, "genbank")
            print(f"Saved → {out_path}")


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
        wb       = openpyxl.load_workbook(path, read_only=True, data_only=True)
        ws       = wb.active
        rows_iter = ws.iter_rows(values_only=True)
        headers  = [str(h) for h in next(rows_iter)]
        rows     = [dict(zip(headers, (str(v) if v is not None else "" for v in row)))
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

    header_font = Font(name="Calibri", bold=True)
    for col_idx, col_name in enumerate(IDT_COLUMNS, start=1):
        cell = ws.cell(row=1, column=col_idx, value=col_name)
        cell.font = header_font

    data_font = Font(name="Calibri")
    for row_idx, row in enumerate(idt_rows, start=2):
        for col_idx, col_name in enumerate(IDT_COLUMNS, start=1):
            cell = ws.cell(row=row_idx, column=col_idx, value=row[col_name])
            cell.font = data_font

    for col_letter, width in IDT_COL_WIDTHS.items():
        ws.column_dimensions[col_letter].width = width

    wb.save(path)
    print(f"Saved → {path}")


def resolve_idt_path(args_output: Path | None, input_path: Path) -> Path:
    """
    Resolve the IDT Excel output path.
      - No --output:             same directory as input, default filename
      - Directory / trailing slash: that directory, default filename
      - Explicit .xlsx path:     used as-is
    """
    default_name = input_path.stem + "_IDT_orderForm.xlsx"
    if args_output is None:
        return input_path.with_name(default_name)
    elif args_output.suffix.lower() == ".xlsx":
        return args_output
    else:
        return args_output / default_name


def resolve_gb_dir(args_gb_output: Path | None, input_path: Path) -> Path:
    """
    Resolve the output directory for .gb files.
    Defaults to the same directory as the input file.
    """
    return args_gb_output if args_gb_output is not None else input_path.parent


# ---------------------------------------------------------------------------
# FASTA helper
# ---------------------------------------------------------------------------

def read_fasta_sequence(path: Path) -> str:
    """
    Read the first sequence from a FASTA file and return it as a plain string.
    Header lines (starting with '>') and whitespace are stripped.

    Parameters
    ----------
    path : Path
        Path to the FASTA file.

    Returns
    -------
    str
        Concatenated sequence of the first FASTA record, whitespace removed.
    """
    seq_lines = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if seq_lines:
                    break       # stop after first record
                continue
            seq_lines.append(line.strip())
    if not seq_lines:
        raise ValueError(f"No sequence found in FASTA file: {path}")
    return "".join(seq_lines)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="sgrna_oligo_designer",
        description=(
            "Generate BsmBI/Golden Gate annealing oligos from an sgRNA spreadsheet, "
            "write a ready-to-upload IDT bulk order Excel file, and optionally produce "
            "annotated full-plasmid GenBank files for sequencing verification."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # --- Input ---
    p.add_argument(
        "input", type=Path,
        help="Input file (.xlsx, .csv, or .tsv) with sgRNA_id and sgRNA_seq columns.",
    )

    # --- IDT Excel output ---
    p.add_argument(
        "--output", "-o", type=Path, default=None,
        help=(
            "Output directory or explicit .xlsx path for the IDT order form. "
            "Default: same directory as input, '<input>_IDT_orderForm.xlsx'."
        ),
    )

    # --- GenBank options ---
    gb = p.add_argument_group("GenBank options (require --vector and biopython)")
    gb.add_argument(
        "--vector", "-v", type=Path, default=None,
        help="Empty vector GenBank file (.gb / .gbk). Required for .gb output.",
    )
    gb.add_argument(
        "--stuffer", default=DEFAULT_STUFFER,
        help=(
            "Placeholder sequence to replace with the sgRNA 20-mer. "
            "Accepts either a raw sequence string or a path to a FASTA file "
            "(auto-detected — if the value is an existing file, it is read as FASTA)."
        ),
    )
    gb.add_argument(
        "--multi-gb", action="store_true",
        help="Write all plasmid records into a single combined .gb file instead of one per sgRNA.",
    )
    gb.add_argument(
        "--gb-output", type=Path, default=None,
        help="Output directory for .gb files. Default: same directory as input.",
    )
    gb.add_argument(
        "--vector-name", default=None,
        help=(
            "Name appended to each .gb output filename: {sgRNA_id}_{vector_name}.gb. "
            "Default: stem of the --vector filename (e.g. 'ipUSEPR' from 'ipUSEPR.gb')."
        ),
    )

    # --- Column / oligo overrides ---
    p.add_argument("--seq-col",    default="sgRNA_seq", help="Column name for sgRNA sequences.")
    p.add_argument("--id-col",     default="sgRNA_id",  help="Column name for sgRNA identifiers.")
    p.add_argument("--top-prefix", default="CACCG",     help="5' prefix for top (sense) oligo.")
    p.add_argument("--bot-prefix", default="AAAC",      help="5' prefix for bottom (antisense) oligo.")
    p.add_argument("--bot-suffix", default="C",         help="3' suffix for bottom (antisense) oligo.")

    return p


def main(argv=None):
    parser = build_parser()
    args   = parser.parse_args(argv)

    input_path: Path = args.input
    if not input_path.exists():
        print(f"Error: input file not found — {input_path}", file=sys.stderr)
        sys.exit(1)

    if args.vector and not args.vector.exists():
        print(f"Error: vector file not found — {args.vector}", file=sys.stderr)
        sys.exit(1)

    # --- Read and process sgRNA table ---
    rows = read_input(input_path)

    # Validate unique sgRNA IDs
    ids        = [row[args.id_col] for row in rows]
    duplicates = [id_ for id_, count in Counter(ids).items() if count > 1]
    if duplicates:
        print(
            f"Error: duplicate sgRNA IDs found: {duplicates}. "
            "Each sgRNA must have a unique identifier.",
            file=sys.stderr,
        )
        sys.exit(1)

    rows = add_cloning_oligos(
        rows,
        seq_col=args.seq_col,
        top_prefix=args.top_prefix,
        bot_prefix=args.bot_prefix,
        bot_suffix=args.bot_suffix,
    )

    # --- IDT Excel (always) ---
    idt_path = resolve_idt_path(args.output, input_path)
    idt_rows = build_idt_rows(rows, id_col=args.id_col)
    write_idt_excel(idt_rows, idt_path)

    # --- GenBank files (only if --vector is provided) ---
    if args.vector:
        # Resolve stuffer — FASTA file or raw string
        stuffer_path = Path(args.stuffer)
        if stuffer_path.exists():
            stuffer_seq = read_fasta_sequence(stuffer_path)
        else:
            stuffer_seq = "".join(args.stuffer.split())

        vector_name   = args.vector_name if args.vector_name else args.vector.stem
        vector_record = read_vector_gb(args.vector)
        records       = build_plasmid_records(
            rows,
            vector_record,
            stuffer_seq=stuffer_seq,
            seq_col=args.seq_col,
            id_col=args.id_col,
        )
        gb_dir = resolve_gb_dir(args.gb_output, input_path)
        write_gb_files(
            records,
            output_dir=gb_dir,
            vector_name=vector_name,
            multi_gb=args.multi_gb,
            combined_stem=input_path.stem + "_plasmids_combined",
        )

    # --- Summary ---
    n = len(rows)
    gb_note = (
        f", {n} plasmid .gb files" if args.vector and not args.multi_gb
        else f", 1 combined plasmid .gb file" if args.vector
        else ""
    )
    print(
        f"Processed {n} sgRNA sequences → "
        f"{n * 2} oligos to order ({n} top + {n} bottom){gb_note}"
    )


if __name__ == "__main__":
    main()

"""
sgrna_oligo_designer.py
------------------------
Core logic for the sgRNA Oligo Designer web app (see app.py): generates annealing
oligo pairs for sgRNA cloning, builds IDT bulk-order rows, and produces annotated
full-plasmid GenBank records for sequencing verification.

This module has no command-line interface; the Shiny app is the only entry point.

Default parameters target BsmBI / U6 vectors such as ipUSEPR:
  - Sense oligo    : 5'-CACCG-{sgRNA_seq}-3'
  - Antisense oligo: 5'-AAAC-{rc(sgRNA_seq)}C-3'
"""

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
    id_col: str  = "sgRNA_id",
    top_prefix: str = "CACCG",
    bot_prefix: str = "AAAC",
    bot_suffix: str = "C",
) -> list[dict]:
    """
    Validate, clean, and add 'oligo_sense' / 'oligo_antisense' keys to each row.

    Each sgRNA sequence is passed through :func:`_validate_sgrna_seq` which:
    - Removes all whitespace (including internal spaces from copy-paste)
    - Uppercases the sequence
    - Validates length (warns if ≠ 20 nt) and character set (ACGT only)
    - Warns if a PAM-like NGG suffix is detected on sequences longer than 20 nt

    The cleaned sequence is written back to ``row[seq_col]`` so downstream
    functions (e.g. GenBank insertion) always use the sanitised version.

    Parameters
    ----------
    rows : list[dict]
        Input rows, each containing at minimum the sgRNA sequence column.
    seq_col : str
        Key holding the raw sgRNA sequences (20-mer, no PAM).
    id_col : str
        Key holding sgRNA identifiers; used in validation warnings / errors.
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
        Same rows with 'oligo_sense' and 'oligo_antisense' added in-place
        and sgRNA sequences cleaned.
    """
    if rows and seq_col not in rows[0]:
        raise ValueError(
            f"Column '{seq_col}' not found. Available: {list(rows[0].keys())}"
        )
    for i, row in enumerate(rows):
        sgrna_id = row.get(id_col, f"row {i + 1}")
        seq      = _validate_sgrna_seq(row[seq_col], sgrna_id)
        row[seq_col]           = seq   # write cleaned sequence back
        row["oligo_sense"]     = top_prefix + seq
        row["oligo_antisense"] = bot_prefix + reverse_complement(seq) + bot_suffix
    return rows


def check_unique_ids(rows: list[dict], id_col: str = "sgRNA_id") -> None:
    """
    Raise ValueError if any sgRNA identifier is missing from the table or repeated.
    """
    if rows and id_col not in rows[0]:
        raise ValueError(
            f"Column '{id_col}' not found. Available: {list(rows[0].keys())}"
        )
    duplicates = [i for i, n in Counter(r[id_col] for r in rows).items() if n > 1]
    if duplicates:
        raise ValueError(
            f"Duplicate sgRNA IDs found: {duplicates}. "
            "Each sgRNA must have a unique identifier."
        )


def build_idt_rows(
    rows: list[dict],
    id_col: str = "sgRNA_id",
) -> list[dict]:
    """
    Reshape oligo rows into IDT bulk upload format.

    Each sgRNA produces two interleaved rows (sense then antisense):
        <sgRNA_id>_s      | oligo_sense     | 25nm | STD
        <sgRNA_id>_a      | oligo_antisense | 25nm | STD

    Parameters
    ----------
    rows : list[dict]
        Rows with sgRNA_id, oligo_sense, and oligo_antisense keys.
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
            "Sequence":     r["oligo_sense"],
            "Scale":        IDT_SCALE,
            "Purification": IDT_PURIF,
        })
        idt_rows.append({
            "Name":         f"{r[id_col]}_a",
            "Sequence":     r["oligo_antisense"],
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


def _validate_sgrna_seq(seq: str, sgrna_id: str) -> str:
    """
    Clean and validate a raw sgRNA sequence.

    Processing order:
      1. Remove ALL whitespace (leading, trailing, and internal) with
         ``"".join(seq.split())`` — handles spaces every 10 bases copied
         from papers, Excel line breaks, and PDF paste artifacts.
      2. Convert to uppercase.
      3. Validate: non-empty, ACGT-only, exactly 20 nt (warn if not).
      4. Warn if length > 20 and sequence ends in NGG — likely PAM
         contamination where the user accidentally included the PAM.

    Parameters
    ----------
    seq : str
        Raw sgRNA sequence from the input table.
    sgrna_id : str
        sgRNA identifier used in error / warning messages.

    Returns
    -------
    str
        Cleaned, uppercase, whitespace-free sequence.

    Raises
    ------
    ValueError
        If the sequence is empty or contains non-ACGT characters.
    """
    seq = "".join(seq.split()).upper()

    if not seq:
        raise ValueError(
            f"sgRNA '{sgrna_id}': sequence is empty. "
            "Check that the sequence column is correct."
        )

    invalid = set(seq) - set("ACGT")
    if invalid:
        raise ValueError(
            f"sgRNA '{sgrna_id}': contains non-ACGT characters: {sorted(invalid)}. "
            "Only standard ACGT bases are supported. Remove ambiguity codes, "
            "spaces, or non-DNA characters before running the script."
        )

    if len(seq) != 20:
        print(
            f"Warning: sgRNA '{sgrna_id}' is {len(seq)} nt "
            f"(expected 20 nt for SpCas9). Proceeding anyway.",
            file=sys.stderr,
        )

    # PAM contamination check — >20 nt ending in NGG
    if len(seq) > 20 and seq[-1] == "G" and seq[-2] == "G":
        print(
            f"Warning: sgRNA '{sgrna_id}' is {len(seq)} nt and ends in "
            f"'{seq[-3:]}' — this looks like an NGG PAM was accidentally "
            "included in the sequence. Trim to 20 nt if so.",
            file=sys.stderr,
        )

    return seq


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

    Features entirely within the stuffer region are removed before shifting —
    this clears any stuffer annotation (misc_feature or similar) that was
    present in the empty vector. A new misc_feature labelled with sgrna_id
    is appended at the insertion site, and all features are sorted by start
    position to follow GenBank convention.

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
        New record with sgRNA inserted, stuffer annotation removed,
        and features sorted by position.
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

    # Remove features entirely within the stuffer region —
    # these annotate the placeholder sequence that has been replaced.
    new_record.features = [
        f for f in new_record.features
        if not (int(f.location.start) >= insert_start
                and int(f.location.end) <= insert_end)
    ]

    # Shift features downstream of the insertion site
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

    # Sort features by start position (GenBank convention)
    new_record.features.sort(key=lambda f: int(f.location.start))

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

    Each record’s LOCUS name (record.name) and accession (record.id) are set to
    '{vector_name}_{sgRNA_id}' so that SnapGene and Benchling display the full
    construct name (e.g. 'ipUSEPR_sg_Luc_474') rather than just the bare sgRNA ID.

    Parameters
    ----------
    records : list[Bio.SeqRecord.SeqRecord]
        Annotated plasmid records to write.
    output_dir : Path
        Directory for output .gb files.
    vector_name : str
        Vector name prepended to each record name and output filename.
        Produces: {vector_name}_{sgRNA_id}.gb
    multi_gb : bool
        If True, write all records into a single combined .gb file.
    combined_stem : str
        Filename stem for the combined file (used when multi_gb=True).
    """
    import copy
    from Bio import SeqIO

    output_dir.mkdir(parents=True, exist_ok=True)

    def _named(record):
        """Return a copy with vector_name prepended to id and name."""
        r             = copy.deepcopy(record)
        full_id       = f"{vector_name}_{record.id}"
        r.id          = full_id
        r.name        = full_id   # LOCUS name shown by SnapGene / Benchling
        r.description = f"{record.id} cloned into {vector_name}"
        return r

    if multi_gb:
        out_path = output_dir / f"{combined_stem}.gb"
        with open(out_path, "w") as fh:
            SeqIO.write([_named(r) for r in records], fh, "genbank")
        print(f"Saved → {out_path}")
    else:
        for record in records:
            named    = _named(record)
            out_path = output_dir / f"{named.id}.gb"
            with open(out_path, "w") as fh:
                SeqIO.write([named], fh, "genbank")
            print(f"Saved → {out_path}")


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def read_input(path: Path) -> list[dict]:
    """
    Read an input spreadsheet or delimited file into a list of row dicts.
    Supports .xlsx / .xlsm / .xls (via openpyxl), .csv, and .tsv.

    Completely empty rows (all values empty or whitespace-only) are filtered
    out — this prevents confusing errors from trailing blank rows that
    openpyxl and csv.DictReader both return for some files.
    """
    suffix = path.suffix.lower()

    def _non_empty(rows: list[dict]) -> list[dict]:
        """Drop rows where every value is empty or whitespace-only."""
        return [
            row for row in rows
            if any(str(v).strip() for v in row.values())
        ]

    if suffix in {".xlsx", ".xlsm", ".xls"}:
        import openpyxl
        wb       = openpyxl.load_workbook(path, read_only=True, data_only=True)
        ws       = wb.active
        rows_iter = ws.iter_rows(values_only=True)
        headers  = [str(h) for h in next(rows_iter)]
        rows     = [dict(zip(headers, (str(v) if v is not None else "" for v in row)))
                    for row in rows_iter]
        wb.close()
        return _non_empty(rows)

    elif suffix == ".csv":
        with open(path, newline="", encoding="utf-8-sig") as fh:
            return _non_empty(list(csv.DictReader(fh)))

    elif suffix == ".tsv":
        with open(path, newline="", encoding="utf-8-sig") as fh:
            return _non_empty(list(csv.DictReader(fh, delimiter="\t")))

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

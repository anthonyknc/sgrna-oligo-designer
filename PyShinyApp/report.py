"""
report.py
---------
Presentation helpers for the web app: result tables, the oligo design template, and the
download bundles. No Shiny code here, so everything is easy to test.
"""

from __future__ import annotations

import csv
import io
import re
import zipfile

from sgrna_oligo_designer import (
    CloningSite,
    build_insert_core,
    design_oligos,
    genbank_string,
    named_record,
    sanitize_name,
    write_idt_excel,
)

SPACER_LEN = 20


def gc_percent(seq: str) -> float:
    """GC content of a DNA sequence in percent (0 for an empty sequence)."""
    return round(100 * sum(c in "GC" for c in seq.upper()) / len(seq), 1) if seq else 0.0


def oligo_table(rows: list[dict], id_col: str = "sgRNA_id", seq_col: str = "sgRNA_seq") -> list[dict]:
    """One display row per sgRNA: spacer, both oligos, their length and GC%, and any notes."""
    out = []
    for r in rows:
        out.append({
            "sgRNA": r[id_col],
            "Spacer (5'→3')": r[seq_col],
            "Sense oligo (5'→3')": r["oligo_sense"],
            "Antisense oligo (5'→3')": r["oligo_antisense"],
            "Length (nt)": f"{len(r['oligo_sense'])} / {len(r['oligo_antisense'])}",
            "GC %": gc_percent(r["oligo_sense"]),
            "Notes": "; ".join(w.split("' ", 1)[-1] if w.startswith("sgRNA '") else w
                               for w in r["warnings"]) or "OK",
        })
    return out


def design_template(site: CloningSite, g_mode: str = "always") -> tuple[str, str]:
    """
    The two oligos for a placeholder spacer of 20 N's, e.g. ``CACCGNNNNNNNNNNNNNNNNNNNN`` and
    ``AAACNNNNNNNNNNNNNNNNNNNNC``. With ``if_needed`` the G is shown, since it is added for most guides.
    """
    core = build_insert_core("N" * SPACER_LEN, "always" if g_mode == "if_needed" else g_mode)
    return design_oligos(core, site)


def split_spacer(template: str) -> tuple[str, str, str]:
    """Split a template into (before, spacer N's, after) for highlighting."""
    m = re.search(r"N{%d}" % SPACER_LEN, template)
    return (template[:m.start()], m.group(0), template[m.end():]) if m else (template, "", "")


def oligos_csv(rows: list[dict], id_col: str = "sgRNA_id", seq_col: str = "sgRNA_seq") -> bytes:
    buf = io.StringIO()
    writer = csv.writer(buf)
    writer.writerow(["sgRNA", "spacer", "sense_oligo", "antisense_oligo", "gc_percent", "notes"])
    for r in oligo_table(rows, id_col, seq_col):
        writer.writerow([r["sgRNA"], r["Spacer (5'→3')"], r["Sense oligo (5'→3')"],
                         r["Antisense oligo (5'→3')"], r["GC %"], r["Notes"]])
    return buf.getvalue().encode()


def _zip(files: dict[str, bytes]) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for name, data in files.items():
            z.writestr(name, data)
    return buf.getvalue()


def genbank_zip(records: list, vector_name: str) -> bytes:
    """A zip with one GenBank file per sgRNA, named '{vector}_{sgRNA}.gb'."""
    files = {}
    for rec in records:
        named = named_record(rec, vector_name)
        files[f"{sanitize_name(named.id)}.gb"] = genbank_string([named]).encode()
    return _zip(files)


def bundle_zip(idt_rows: list[dict], rows: list[dict], records: list, vector_name: str,
               combined_stem: str, tmp_xlsx_path, id_col: str = "sgRNA_id",
               seq_col: str = "sgRNA_seq") -> bytes:
    """Everything in one zip: IDT order form, oligo table (CSV), combined GenBank, per-sgRNA GenBank."""
    write_idt_excel(idt_rows, tmp_xlsx_path)
    named = [named_record(r, vector_name) for r in records]
    files = {
        "IDT_orderForm.xlsx": tmp_xlsx_path.read_bytes(),
        "oligos.csv": oligos_csv(rows, id_col, seq_col),
        f"{sanitize_name(combined_stem, 'combined_plasmids')}.gb": genbank_string(named).encode(),
    }
    for rec in named:
        files[f"plasmids/{sanitize_name(rec.id)}.gb"] = genbank_string([rec]).encode()
    return _zip(files)

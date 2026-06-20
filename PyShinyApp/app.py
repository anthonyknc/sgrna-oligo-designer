"""
app.py  —  sgRNA Oligo Designer v1.0.0
---------------------------------------
Core shiny mode. Run: shiny run app.py --reload
Deploy:        rsconnect deploy shiny . --name <account> --title sgrna-oligo-designer
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pandas as pd
import shinyswatch
from shiny import App, reactive, render, req
from shiny import ui as sui

from sgrna_oligo_designer import (
    add_cloning_oligos,
    build_idt_rows,
    build_plasmid_records,
    read_fasta_sequence,
    read_input,
    read_vector_gb,
    write_gb_files,
    write_idt_excel,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_VERSION = "v1.0.0"
_GITHUB_URL = "https://github.com/anthonyknc/sgrna-oligo-designer"

_GITHUB_SVG = """
<svg xmlns="http://www.w3.org/2000/svg" width="22" height="22"
     viewBox="0 0 24 24" fill="currentColor" aria-label="GitHub">
  <path d="M12 0c-6.626 0-12 5.373-12 12 0 5.302 3.438 9.8 8.207
  11.387.599.111.793-.261.793-.577v-2.234c-3.338.726-4.033-1.416-4.033
  -1.416-.546-1.387-1.333-1.756-1.333-1.756-1.089-.745.083-.729.083-.729
  1.205.084 1.839 1.237 1.839 1.237 1.07 1.834 2.807 1.304 3.492.997
  .107-.775.418-1.305.762-1.604-2.665-.305-5.467-1.334-5.467-5.931
  0-1.311.469-2.381 1.236-3.221-.124-.303-.535-1.524.117-3.176 0 0
  1.008-.322 3.301 1.23.957-.266 1.983-.399 3.003-.404 1.02.005 2.047
  .138 3.006.404 2.291-1.552 3.297-1.23 3.297-1.23.653 1.653.242 2.874
  .118 3.176.77.84 1.235 1.911 1.235 3.221 0 4.609-2.807 5.624-5.479
  5.921.43.372.823 1.102.823 2.222v3.293c0 .319.192.694.801.576
  4.765-1.589 8.199-6.086 8.199-11.386 0-6.627-5.373-12-12-12z"/>
</svg>
"""

# Bootstrap 5.3 native dark mode toggle via data-bs-theme attribute
_DARK_MODE_JS = """
Shiny.addCustomMessageHandler('set-theme', function(theme) {
    document.documentElement.setAttribute('data-bs-theme', theme);
});
"""

# ---------------------------------------------------------------------------
# Title bar: name + version + GitHub icon
# ---------------------------------------------------------------------------

_TITLE = sui.tags.div(
    sui.tags.span("sgRNA Oligo Designer", style="font-weight:600;"),
    sui.tags.span(
        _VERSION,
        style="font-size:0.7em; color:#6c757d; margin-left:0.5em; vertical-align:middle;",
    ),
    sui.tags.a(
        sui.HTML(_GITHUB_SVG),
        href=_GITHUB_URL,
        target="_blank",
        title="View on GitHub",
        style="margin-left:0.75em; color:inherit; vertical-align:middle; opacity:0.75;",
    ),
    style="display:flex; align-items:center;",
)

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

app_ui = sui.page_sidebar(

    sui.sidebar(

        sui.h6("Step 1 — sgRNA table"),
        sui.input_file(
            "sgrna_file",
            "Upload sgRNA table (.xlsx, .csv, or .tsv)",
            accept=[".xlsx", ".xlsm", ".csv", ".tsv"],
        ),

        sui.hr(),
        sui.h6("Step 2 — GenBank output (optional)"),
        sui.input_file(
            "vector_file",
            "Upload empty vector (.gb / .gbk)",
            accept=[".gb", ".gbk"],
        ),
        sui.input_text(
            "stuffer_text",
            "Stuffer sequence",
            placeholder="e.g. GAGACGTGACCGTCTCT",
        ),
        sui.input_file(
            "stuffer_file",
            "Or upload stuffer FASTA (.fa)",
            accept=[".fa", ".fasta"],
        ),
        # Stuffer source indicator — updates reactively without Run
        sui.output_ui("stuffer_source"),
        sui.input_text(
            "vector_name",
            "Vector name (optional)",
            placeholder="Default: vector filename stem",
        ),
        sui.input_text(
            "combined_name",
            "Combined .gb filename stem",
            value="combined_plasmids",
        ),

        sui.hr(),
        sui.h6("Advanced options"),
        sui.input_text("seq_col",          "sgRNA sequence column", value="sgRNA_seq"),
        sui.input_text("id_col",           "sgRNA ID column",        value="sgRNA_id"),
        sui.input_text("sense_prefix",     "Sense oligo prefix",     value="CACCG"),
        sui.input_text("antisense_prefix", "Antisense oligo prefix", value="AAAC"),
        sui.input_text("antisense_suffix", "Antisense oligo suffix", value="C"),

        sui.hr(),
        sui.input_action_button("run_btn", "Run", class_="btn-primary w-100"),

        sui.hr(),
        sui.input_switch("dark_mode", "🌙  Dark mode", value=False),

        width=320,
    ),

    # Dark mode JS handler
    sui.tags.script(_DARK_MODE_JS),

    # --- Main panel ---
    sui.h4("Summary"),
    sui.output_text("summary"),

    sui.hr(),
    sui.h4("IDT order form preview"),
    sui.output_table("idt_table"),

    sui.hr(),
    sui.h4("Downloads"),
    sui.h6("IDT Bulk Order Form (.xlsx)"),
    sui.download_button("download_idt", "Download IDT Order Form"),
    sui.h6("Combined GenBank File (.gb)"),
    sui.download_button("download_gb", "Download Combined GenBank"),

    title=_TITLE,
    theme=shinyswatch.theme.flatly(),
)


# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

def server(input, output, session):

    # --- Dark mode toggle ---

    @reactive.effect
    @reactive.event(input.dark_mode)
    async def apply_theme():
        theme = "dark" if input.dark_mode() else "light"
        await session.send_custom_message("set-theme", theme)

    # --- Stuffer source indicator (reactive, no Run needed) ---

    @output
    @render.ui
    def stuffer_source():
        if input.stuffer_file():
            return sui.tags.small(
                "✓ Using uploaded FASTA file",
                style="color:#2ecc71; font-weight:500;",
            )
        elif input.stuffer_text().strip():
            return sui.tags.small(
                "✓ Using manually entered sequence",
                style="color:#2ecc71; font-weight:500;",
            )
        else:
            return sui.tags.small(
                "⚠ No stuffer — GenBank output disabled",
                style="color:#f39c12; font-weight:500;",
            )

    # --- Reactive calculations ---

    @reactive.calc
    @reactive.event(input.run_btn)
    def processed_rows() -> list[dict]:
        req(input.sgrna_file())
        file_info = input.sgrna_file()[0]
        rows = read_input(Path(file_info["datapath"]))
        return add_cloning_oligos(
            rows,
            seq_col=input.seq_col(),
            top_prefix=input.sense_prefix(),
            bot_prefix=input.antisense_prefix(),
            bot_suffix=input.antisense_suffix(),
        )

    @reactive.calc
    @reactive.event(input.run_btn)
    def idt_rows() -> list[dict]:
        return build_idt_rows(processed_rows(), id_col=input.id_col())

    @reactive.calc
    @reactive.event(input.run_btn)
    def plasmid_records():
        req(input.vector_file())
        vector_info = input.vector_file()[0]
        vector_record = read_vector_gb(Path(vector_info["datapath"]))
        if input.stuffer_file():
            stuffer_seq = read_fasta_sequence(Path(input.stuffer_file()[0]["datapath"]))
        else:
            stuffer_seq = "".join(input.stuffer_text().split())
        req(stuffer_seq)
        return build_plasmid_records(
            processed_rows(),
            vector_record,
            stuffer_seq=stuffer_seq,
            seq_col=input.seq_col(),
            id_col=input.id_col(),
        )

    # --- Outputs ---

    @output
    @render.text
    def summary():
        req(processed_rows())
        n = len(processed_rows())
        try:
            gb_note = " | GenBank file ready" if plasmid_records() else ""
        except Exception:
            gb_note = ""
        return (
            f"Processed {n} sgRNA sequences "
            f"→ {n * 2} oligos ({n} sense + {n} antisense)"
            f"{gb_note}"
        )

    @output
    @render.table
    def idt_table():
        req(idt_rows())
        return pd.DataFrame(idt_rows())

    # --- Downloads ---

    @render.download(filename="IDT_orderForm.xlsx")
    def download_idt():
        req(idt_rows())
        with tempfile.NamedTemporaryFile(suffix=".xlsx", delete=False) as tmp:
            tmp_path = Path(tmp.name)
        try:
            write_idt_excel(idt_rows(), tmp_path)
            yield tmp_path.read_bytes()
        finally:
            tmp_path.unlink(missing_ok=True)

    @render.download(filename="combined_plasmids.gb")
    def download_gb():
        req(plasmid_records())
        vector_name = (
            input.vector_name()
            or Path(input.vector_file()[0]["name"]).stem
        )
        combined_stem = input.combined_name() or "combined_plasmids"
        with tempfile.TemporaryDirectory() as tmp_dir:
            write_gb_files(
                plasmid_records(),
                output_dir=Path(tmp_dir),
                vector_name=vector_name,
                multi_gb=True,
                combined_stem=combined_stem,
            )
            gb_file = Path(tmp_dir) / f"{combined_stem}.gb"
            yield gb_file.read_bytes()


# ---------------------------------------------------------------------------
# App
# ---------------------------------------------------------------------------

app = App(app_ui, server)

"""
app.py  —  sgRNA Oligo Designer (core shiny mode)
--------------------------------------------------
Core shiny mode is used here for maximum compatibility.
The reactive model is identical to R/Shiny: app_ui = ui.R, server = server.R.

Run locally:
    shiny run app.py --reload

Deploy to shinyapps.io:
    rsconnect deploy shiny . --name <your-account> --title sgrna-oligo-designer
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pandas as pd
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
# UI  (equivalent to ui.R in R/Shiny)
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
        width=320,
    ),

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

    title="sgRNA Oligo Designer",
)


# ---------------------------------------------------------------------------
# Server  (equivalent to server.R in R/Shiny)
# ---------------------------------------------------------------------------

def server(input, output, session):

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

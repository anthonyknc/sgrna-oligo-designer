"""
app.py  —  sgRNA Oligo Designer
--------------------------------
Shiny for Python app (Core mode). Run: shiny run app.py --reload
All design logic lives in sgrna_oligo_designer.py; presentation helpers are in report.py.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pandas as pd
import shinyswatch
from faicons import icon_svg
from shiny import App, reactive, render, req, ui
from shiny.types import SilentException

from report import (
    TABLE_COLUMNS,
    bundle_zip,
    design_template,
    genbank_zip,
    oligo_table,
    oligos_csv,
    split_spacer,
)
from sgrna_oligo_designer import (
    ENZYMES,
    add_cloning_oligos,
    all_enzymes,
    build_idt_rows,
    build_plasmid_records,
    find_cloning_site,
    genbank_string,
    named_record,
    read_input,
    read_vector_gb,
    sanitize_name,
    single_cutters,
    suggest_cloning_sites,
    validate_ids,
    write_idt_excel,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_VERSION      = "v2.1.0"
_LAST_UPDATE  = "2026-11-23"   # YYYY-MM-DD; bump together with _VERSION on each release
_GITHUB_URL   = "https://github.com/anthonyknc/sgrna-oligo-designer"
_HERE         = Path(__file__).parent
_EXAMPLE_XLSX = _HERE / "example_inputs" / "example_sgRNA_seqs.xlsx"

# Vectors bundled with the app: key -> (label, GenBank file)
_VECTORS = {
    "ipUSEPR":          ("ipUSEPR",          _HERE / "vectors" / "ipUSEPR.gb"),
    "lentiCRISPRv2":    ("lentiCRISPR v2",   _HERE / "vectors" / "lentiCRISPRv2.gb"),
    "LentiCRISPRv2GFP": ("LentiCRISPRv2GFP", _HERE / "vectors" / "LentiCRISPRv2GFP.gb"),
    "lentiGuide-Puro":  ("lentiGuide-Puro (Addgene 52963)", _HERE / "vectors" / "lentiGuide-Puro.gb"),
    "pX330":            ("pX330 (Addgene 42230)", _HERE / "vectors" / "pX330.gb"),
    "pU6-BbsI-CBh-Cas9-T2A-mCherry-H1":
                        ("pU6-(BbsI)_CBh-Cas9-T2A-mCherry-H1 (Addgene 64217)",
                         _HERE / "vectors" / "pU6-BbsI_CBh-Cas9-T2A-mCherry-H1.gb"),
    "phH1-gRNA":        ("phH1-gRNA, H1 promoter (Addgene 53186)", _HERE / "vectors" / "phH1-gRNA.gb"),
    "pU6-sgRNA-EF1a-puro-T2A-BFP":
                        ("pU6-sgRNA EF1Alpha-puro-T2A-BFP, BstXI + BlpI (Addgene 60955)",
                         _HERE / "vectors" / "pU6-sgRNA_EF1Alpha-puro-T2A-BFP.gb"),
    "pMCB306":          ("pMCB306, mouse U6, BstXI + BlpI (Addgene 89360)", _HERE / "vectors" / "pMCB306.gb"),
}
_UPLOAD = "upload"
_OTHER = "__other__"      # enzyme dropdown value that reveals the full enzyme search

# Same-site enzymes are shown together; either name works (Esp3I and BpiI are also under "Other").
_ENZYME_LABELS = {"BsmBI": "BsmBI / Esp3I", "BbsI": "BbsI / BpiI"}

_G_CHOICES = {
    "always":    "Always add G (standard)",
    "if_needed": "Only if the sgRNA doesn't start with G",
    "never":     "Never (non-U6 promoter)",
}
_ORIENT_CHOICES = {
    "auto":    "Auto-detect",
    "forward": "Forward (as in the GenBank file)",
    "reverse": "Reverse",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _info(title: str, *content, placement: str = "right"):
    """
    A small '?' button that opens a popover when clicked (click only, so it also works on
    touch screens; click again or click elsewhere to close). A real <button>, so it is
    reachable by keyboard and does not trigger the surrounding label's input.
    """
    return ui.popover(
        ui.tags.button("?", type="button", class_="info-btn",
                        **{"aria-label": f"More information: {title}"}),
        *content,
        title=title,
        placement=placement,
        options={"trigger": "click"},
    )


def _label(text: str, *tip, title: str | None = None):
    """Input label followed by a click-to-open '?' popover."""
    return ui.tags.span(text, _info(title or text, *tip))


def _link(text: str, href: str):
    return ui.tags.a(text, href=href, target="_blank", rel="noopener")


_CRISPICK = "https://portals.broadinstitute.org/gppx/crispick/public"

_TIP_TABLE = (
    "A spreadsheet (.xlsx, .csv or .tsv) with one row per sgRNA and two columns: ",
    ui.tags.code("sgRNA_id"), " (unique, no spaces) and ", ui.tags.code("sgRNA_seq"),
    " (the 20-nt target sequence, without the PAM or any overhangs). Extra columns are ignored, "
    "and spaces and lowercase are cleaned automatically. Design guides with ",
    _link("CRISPick", _CRISPICK), ". More under Guide → Getting started.",
)
_TIP_VECTOR = (
    "The empty vector you will clone into. Nine popular vectors are built in, and the enzyme(s) "
    "that suit each one are chosen for you. Choose “Upload my own” for another circular GenBank file.",
)
_TIP_ENZYME = (
    "The Type IIS enzyme used for sgRNA cloning. The app cuts the vector in silico at the "
    "enzyme's two sites and reads the overhangs from the cuts, so you never type them. The vector "
    "must have exactly two sites for the enzyme you pick. Choose “Other…” to search every "
    "commercially available Type IIS enzyme with a 3–5 nt overhang (e.g. BspQI, BtgZI, FokI). "
    "For an enzyme not on the common list, check that it suits your vector (for example "
    "methylation sensitivity) and follow the supplier's protocol.",
)
_TIP_MODE = (
    "Most sgRNA vectors are cut by one Type IIS enzyme at two sites. Some, such as BstXI + BlpI or "
    "AgeI + EcoRI, are cut by two different enzymes instead; choose “Two different enzymes” for those "
    "and pick them from the enzymes that cut your vector exactly once.",
)
_TIP_FLANK = (
    "The constant bases next to the spacer in your vector's oligo design (for BstXI/BlpI vectors: "
    "TTG before it and GTTTAAGAGC after it). Leave blank to use the values worked out from the "
    "vector; type bases to override them. Check them against your vector's cloning protocol.",
)
_TIP_G = (
    "The human U6 promoter starts transcription best at a G. ‘Always’ adds a G to every sgRNA "
    "(CACCG-sgRNA), the standard in the lentiCRISPR and Broad protocols, even when the sgRNA "
    "already starts with G. See ",
    _link("Ran et al., 2013", "https://pubmed.ncbi.nlm.nih.gov/24157548/"),
    ". Details under Guide → How it works.",
)
_TIP_SITE = (
    "The 4-nt overhangs your oligos must carry, read from the cut ends of the vector, and how many "
    "bases between the cuts are replaced by the sgRNA. Orientation is the strand of the GenBank "
    "file the sgRNA is transcribed from.",
)


def _guide(name: str):
    """Render one of the help/*.md files (the single source of truth for the guide)."""
    return ui.markdown((_HERE / "help" / f"{name}.md").read_text(encoding="utf-8"))


# Close any open info popover on an outside click or Escape. (Bootstrap's click trigger only
# toggles; its 'focus' trigger would break in Safari, which doesn't focus buttons on click.)
_INFO_JS = """
(function () {
  function closeOpen(except, refocus) {
    document.querySelectorAll('.popover.show').forEach(function (pop) {
      var trigger = document.querySelector('[aria-describedby="' + pop.id + '"]');
      if (!trigger || trigger === except) return;
      var inst = window.bootstrap && bootstrap.Popover.getInstance(trigger);
      if (inst) inst.hide();
      if (refocus) trigger.focus();
    });
  }
  document.addEventListener('click', function (e) {
    if (e.target.closest('.popover')) return;
    closeOpen(e.target.closest('.info-btn'));
  });
  document.addEventListener('keydown', function (e) { if (e.key === 'Escape') closeOpen(null, true); });
})();
"""


# ---------------------------------------------------------------------------
# Styling (uses Bootstrap variables, so it follows light/dark mode automatically)
# ---------------------------------------------------------------------------

def _fa(name: str, **kw):
    """Font Awesome icon that sits nicely inline with text."""
    return icon_svg(name, margin_left="0", margin_right="0.45em", **kw)


_CSS = """
:root { --brand: #0f766e; --brand-bg: rgba(15,118,110,.12); }
[data-bs-theme="dark"] { --brand: #5eead4; --brand-bg: rgba(94,234,212,.16); }
body { font-feature-settings: "cv11", "ss01"; }
.navbar-brand { font-weight: 650; letter-spacing: -0.01em; display: inline-flex; align-items: center; }
.navbar-brand .badge { font-size: 0.62rem; font-weight: 600; margin-left: 0.6rem; vertical-align: middle; }
.card { border: 1px solid var(--bs-border-color); border-radius: 0.85rem;
        box-shadow: 0 1px 2px rgba(0,0,0,.04); }
.card-header { font-weight: 600; display: flex; align-items: center; gap: .5rem;
               background: transparent; border-bottom: 1px solid var(--bs-border-color); }
.card-header .header-actions { margin-left: auto; display: flex; flex-wrap: wrap; gap: .4rem; }
.bslib-sidebar-layout > .sidebar { background: var(--bs-tertiary-bg); }
.sidebar .accordion-button { font-weight: 600; padding: .7rem .9rem; }
.sidebar .accordion-body { padding: .8rem .9rem 1rem; }
.sidebar .form-label, .sidebar label { font-size: .86rem; }
#run_btn { font-weight: 600; letter-spacing: .01em; }
.shiny-busy #run_btn { opacity: .65; cursor: wait; }
.mono, code.mono, .tpl code { font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; }
.chip { display: inline-flex; align-items: center; gap: .35rem; padding: .2rem .65rem; margin: .15rem .3rem .15rem 0;
        border-radius: 999px; font-size: .82rem; background: var(--bs-tertiary-bg);
        border: 1px solid var(--bs-border-color); }
.chip .k { color: var(--bs-secondary-color); }
.chip code { color: inherit; background: transparent; padding: 0; font-weight: 600; }
.status-ok { border-left: 4px solid var(--bs-success); }
.status-bad { border-left: 4px solid var(--bs-danger); }
.tpl { margin-top: .6rem; padding: .7rem .9rem; border-radius: .6rem; background: var(--bs-tertiary-bg);
       border: 1px dashed var(--bs-border-color); }
.tpl-row { display: flex; gap: .6rem; align-items: baseline; flex-wrap: wrap; }
.tpl-label { min-width: 5.2rem; font-size: .78rem; color: var(--bs-secondary-color); text-transform: uppercase;
             letter-spacing: .04em; }
.tpl code { font-size: .92rem; color: var(--bs-body-color); background: transparent; padding: 0; word-break: break-all; }
.tpl-spacer { color: var(--brand); font-weight: 700; background: var(--brand-bg); border-radius: .25rem; padding: 0 .15rem; }
.empty-state { text-align: center; padding: 3rem 1rem; color: var(--bs-secondary-color); }
.empty-state svg { width: 2.6rem; height: 2.6rem; color: var(--brand); margin-bottom: .8rem; }
.bslib-value-box { min-height: 0 !important; }
.bslib-value-box .value-box-showcase svg { width: 2rem; height: 2rem; }
.bslib-value-box .value-box-title { font-size: .8rem; text-transform: uppercase; letter-spacing: .05em; opacity: .9; }
.bslib-value-box .value-box-value { font-size: 1.7rem; font-weight: 700; }
.empty-state h5 { color: var(--bs-body-color); }
.empty-state ol { display: inline-block; text-align: left; margin-top: .4rem; }
.info-btn { display: inline-flex; align-items: center; justify-content: center; width: 1.5rem; height: 1.5rem;
            margin-left: 0.4rem; padding: 0; border: 1px solid var(--bs-secondary-color); border-radius: 50%;
            background: transparent; color: var(--bs-secondary-color); font-size: 0.7rem; font-weight: 700;
            line-height: 1; cursor: pointer; vertical-align: middle; }
.info-btn:hover, .info-btn:focus-visible { color: var(--bs-body-bg); background: var(--bs-secondary-color); }
.popover { max-width: 320px; }
.guide { max-width: 54rem; }
.guide h1 { font-size: 1.7rem; margin-bottom: 1rem; }
.guide h2 { font-size: 1.25rem; margin-top: 1.9rem; }
.guide table { width: 100%; margin: 1rem 0; border-collapse: collapse; }
.guide th, .guide td { padding: .45rem .7rem; border-bottom: 1px solid var(--bs-border-color); vertical-align: top; }
.guide th { font-weight: 600; text-align: left; background: var(--bs-tertiary-bg); text-transform: none; }
.guide .nav-pills .nav-link, .nav-pills .nav-link:not(.active) { color: var(--bs-body-color); }
.guide code, .guide pre { overflow-wrap: anywhere; }
.hint { color: var(--bs-secondary-color); }
.app-footer { text-align: center; color: var(--bs-secondary-color); font-size: .8rem; padding: 1.5rem 0 1rem; }
"""


# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

_TITLE = ui.tags.span(
    _fa("dna"), "sgRNA Oligo Designer",
    ui.tags.span(_VERSION, class_="badge text-bg-secondary rounded-pill"),
)


def _sidebar():
    return ui.sidebar(
        ui.input_action_button("example_btn", [_fa("flask"), "Load example sgRNAs"],
                               class_="btn-outline-primary w-100"),
        ui.accordion(
            ui.accordion_panel(
                "1 · sgRNA table",
                ui.input_file(
                    "sgrna_file",
                    _label("Upload sgRNA table", *_TIP_TABLE, title="sgRNA table"),
                    accept=[".xlsx", ".xlsm", ".csv", ".tsv"],
                    button_label="Browse…", placeholder=".xlsx, .csv or .tsv",
                ),
                value="table", icon=_fa("table"),
            ),
            ui.accordion_panel(
                "2 · Vector & enzymes",
                ui.input_select(
                    "vector_choice",
                    _label("Vector", *_TIP_VECTOR),
                    {**{k: v[0] for k, v in _VECTORS.items()}, _UPLOAD: "Upload my own vector (.gb)…"},
                ),
                ui.panel_conditional(
                    f"input.vector_choice == '{_UPLOAD}'",
                    ui.input_file("vector_file", "Upload empty vector (.gb / .gbk)",
                                  accept=[".gb", ".gbk"], button_label="Browse…", placeholder=".gb or .gbk"),
                ),
                ui.input_radio_buttons(
                    "enzyme_mode",
                    _label("Cloning enzymes", *_TIP_MODE),
                    {"one": "One enzyme (cuts the vector twice)", "two": "Two different enzymes (each cuts once)"},
                    selected="one",
                ),
                ui.panel_conditional(
                    "input.enzyme_mode == 'one'",
                    ui.input_select(
                        "enzyme",
                        _label("Restriction enzyme", *_TIP_ENZYME),
                        {**{e: _ENZYME_LABELS.get(e, e) for e in ENZYMES}, _OTHER: "Other…"},
                    ),
                    ui.panel_conditional(
                        f"input.enzyme == '{_OTHER}'",
                        ui.input_selectize(
                            "other_enzyme", "Search all enzymes",
                            {n: f"{n}  ({e.site}, {e.overhang_len}-nt overhang)" for n, e in all_enzymes().items()},
                            selected="SapI",
                        ),
                    ),
                ),
                ui.panel_conditional(
                    "input.enzyme_mode == 'two'",
                    ui.input_selectize("enzyme_a", "Enzyme 1 (cuts once)", {"": "Choose an enzyme…"}),
                    ui.input_selectize("enzyme_b", "Enzyme 2 (cuts once)", {"": "Choose an enzyme…"}),
                ),
                value="vector", icon=_fa("scissors"),
            ),
            ui.accordion_panel(
                "3 · Advanced options",
                ui.input_text(
                    "id_col",
                    _label("sgRNA ID column",
                           "Column name in your input file containing unique sgRNA identifiers."),
                    value="sgRNA_id",
                ),
                ui.input_text(
                    "seq_col",
                    _label("sgRNA sequence column",
                           "Column name in your input file containing the 20-nt sgRNA target sequences."),
                    value="sgRNA_seq",
                ),
                ui.panel_conditional(
                    "input.enzyme_mode == 'two'",
                    ui.input_text("left_extra", _label("Bases before the spacer", *_TIP_FLANK, title="Flanking bases"),
                                  placeholder="auto"),
                    ui.input_text("right_extra", _label("Bases after the spacer", *_TIP_FLANK, title="Flanking bases"),
                                  placeholder="auto"),
                ),
                ui.input_radio_buttons("g_mode", _label("5' G for U6 transcription", *_TIP_G),
                                       _G_CHOICES, selected="always"),
                ui.input_select(
                    "orientation",
                    _label("sgRNA orientation in the vector",
                           "Which strand of the GenBank file the sgRNA is transcribed from. "
                           "Auto-detect looks for the start of the Cas9 scaffold (GTTT)."),
                    _ORIENT_CHOICES,
                ),
                ui.input_text(
                    "vector_name",
                    _label("Vector name (optional)",
                           "Name prepended to each plasmid name in the GenBank files. "
                           "Defaults to the vector name or the uploaded file's name."),
                    placeholder="Default: vector name",
                ),
                ui.input_text(
                    "combined_name",
                    _label("Combined .gb filename stem", "The combined GenBank file will be named <stem>.gb."),
                    value="combined_plasmids",
                ),
                value="advanced", icon=_fa("gear"),
            ),
            id="sidebar_acc", open=["table", "vector"], multiple=True,
        ),
        ui.input_action_button("run_btn", [_fa("play"), "Design oligos"], class_="btn-primary btn-lg w-100"),
        width=390, gap="1rem",
    )


def _design_page():
    return ui.layout_sidebar(
        _sidebar(),
        ui.card(
            ui.card_header(_fa("scissors"), "Cloning site", _info("Cloning site", *_TIP_SITE, placement="bottom")),
            ui.output_ui("site_info"),
        ),
        ui.output_ui("stale_notice"),
        ui.output_ui("results_ui"),
        fillable=False, fill=False,
    )


def _guide_page():
    def page(name: str, title: str):
        return ui.nav_panel(title, ui.div(_guide(name), class_="guide"))
    return ui.div(
        ui.navset_pill_list(
            page("getting_started", "Getting started"),
            page("how_it_works", "How it works"),
            page("troubleshooting", "Troubleshooting"),
            page("citing", "Citing"),
            widths=(3, 9),
        ),
        style="padding-top: 1rem;",
    )


app_ui = ui.page_navbar(
    ui.nav_panel("Design", _design_page(), icon=_fa("dna"), value="design"),
    ui.nav_panel("Guide", _guide_page(), icon=_fa("book-open"), value="guide"),
    ui.nav_spacer(),
    ui.nav_control(ui.tags.a(icon_svg("github", margin_left="0", margin_right="0"), href=_GITHUB_URL,
                             target="_blank", rel="noopener", class_="nav-link", title="View on GitHub",
                             **{"aria-label": "View on GitHub"})),
    ui.nav_control(ui.input_dark_mode(id="dark_mode")),
    title=_TITLE,
    id="main_nav",
    footer=ui.div(
        ui.tags.hr(),
        ui.p("© 2026 Anthony Chan · ", ui.tags.a("GitHub", href=_GITHUB_URL, target="_blank", rel="noopener"),
             " · MIT License"),
        ui.p(f"{_VERSION} · Last Update: {_LAST_UPDATE}", class_="mb-0"),
        class_="app-footer",
    ),
    header=ui.head_content(ui.tags.style(_CSS), ui.tags.script(_INFO_JS)),
    theme=shinyswatch.theme.zephyr(),
    window_title="sgRNA Oligo Designer",
    lang="en",
    fillable=False,
)


# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

def _chip(label: str, value: str):
    return ui.span(ui.span(label, class_="k"), ui.tags.code(value), class_="chip")


def _empty_state():
    return ui.card(
        ui.div(
            icon_svg("dna", margin_left="0", margin_right="0"),
            ui.h5("Design your oligos"),
            ui.tags.ol(
                ui.tags.li("Upload your sgRNA table, or load the example."),
                ui.tags.li("Choose your vector and cloning enzymes."),
                ui.tags.li(["Click ", ui.tags.strong("Design oligos"), "."]),
            ),
            class_="empty-state",
        )
    )


def _alert(kind: str, *children):
    return ui.div(*children, class_=f"alert alert-{kind}", role="alert")


def _error_alert(e: Exception):
    """Our own ValueErrors are written for the user; anything else gets a generic wrapper."""
    if isinstance(e, ValueError):
        return _alert("danger", _fa("circle-xmark"), str(e))
    return _alert("danger", _fa("circle-xmark"),
                  f"Something unexpected went wrong ({type(e).__name__}: {e}). Please check your files, "
                  "and report the problem on GitHub if it persists.")


def server(input, output, session):

    # Example mode flag
    _use_example = reactive.value(False)

    @reactive.effect
    @reactive.event(input.example_btn)
    def _set_example():
        _use_example.set(True)

    # --- Stale results: any setting change after a run invalidates the shown results ---

    _ran   = reactive.value(False)
    _stale = reactive.value(False)

    @reactive.effect
    @reactive.event(input.sgrna_file, ignore_init=True)
    def _uploaded_table():
        _use_example.set(False)      # a newly uploaded table replaces the example

    @reactive.effect
    @reactive.event(input.run_btn, input.example_btn, ignore_init=True)
    def _mark_fresh():
        _ran.set(True)
        _stale.set(False)

    @reactive.effect
    @reactive.event(
        input.vector_choice, input.vector_file, input.enzyme, input.other_enzyme, input.orientation,
        input.enzyme_mode, input.enzyme_a, input.enzyme_b, input.left_extra, input.right_extra,
        input.g_mode, input.seq_col, input.id_col, input.sgrna_file,
        ignore_init=True,
    )
    def _mark_stale():
        if _ran():
            _stale.set(True)

    @render.ui
    def stale_notice():
        if _stale():
            return _alert("info", _fa("circle-info"),
                          "Settings changed since the last run. Click ", ui.tags.strong("Design oligos"),
                          " to update the results.")

    # --- Vector and cloning site (live: updates as soon as the selection changes) ---

    @reactive.calc
    def vector_record():
        if input.vector_choice() == _UPLOAD:
            req(input.vector_file())
            return read_vector_gb(Path(input.vector_file()[0]["datapath"]))
        return read_vector_gb(_VECTORS[input.vector_choice()][1])

    @reactive.calc
    def vector_label() -> str:
        if input.vector_name().strip():
            return sanitize_name(input.vector_name(), "vector")
        if input.vector_choice() == _UPLOAD:
            return sanitize_name(Path(input.vector_file()[0]["name"]).stem, "vector")
        return input.vector_choice()

    @reactive.calc
    def enzyme_name() -> str:
        return input.other_enzyme() if input.enzyme() == _OTHER else input.enzyme()

    def _two_enzymes():
        """The two chosen enzymes, or None until both are picked."""
        a, b = input.enzyme_a(), input.enzyme_b()
        return (a, b) if a and b else None

    def _missing_setup() -> str | None:
        """What is missing from the vector/enzyme settings, in plain words."""
        if input.vector_choice() == _UPLOAD and not input.vector_file():
            return "Upload a vector GenBank file, or choose one of the built-in vectors."
        if input.enzyme_mode() == "two" and _two_enzymes() is None:
            return "Choose the two enzymes."
        return None

    def _missing_input() -> str | None:
        """What the user still has to provide before oligos can be designed."""
        if not _use_example() and not input.sgrna_file():
            return "Upload an sgRNA table, or click Load example sgRNAs."
        return _missing_setup()

    _applied = {"vec": None}      # the vector the suggested enzymes were last applied for

    @reactive.effect
    def _auto_pick():
        """When a vector is chosen, select the enzyme(s) that cut next to its Pol III promoter."""
        found = suggestions()
        if not found:
            return
        first = found[0]
        if first.enzyme_b:
            ui.update_radio_buttons("enzyme_mode", selected="two")
        elif first.enzyme in ENZYMES:
            ui.update_radio_buttons("enzyme_mode", selected="one")
            ui.update_select("enzyme", selected=first.enzyme)
        else:
            ui.update_radio_buttons("enzyme_mode", selected="one")
            ui.update_select("enzyme", selected=_OTHER)
            ui.update_selectize("other_enzyme", selected=first.enzyme)

    @reactive.effect
    def _offer_single_cutters():
        """In two-enzyme mode, offer only the enzymes that cut this vector exactly once."""
        if input.enzyme_mode() != "two":
            return
        try:
            record = vector_record()
            cutters = single_cutters(record)
        except Exception:
            return
        with reactive.isolate():
            keep = {"enzyme_a": input.enzyme_a(), "enzyme_b": input.enzyme_b()}
            suggested = next((x for x in suggestions() if x.enzyme_b), None)
        if _applied["vec"] != id(record):             # a new vector: use its suggested pair, not the old one
            _applied["vec"] = id(record)
            keep = {"enzyme_a": suggested.enzyme if suggested else "",
                    "enzyme_b": suggested.enzyme_b if suggested else ""}
        for id_, current in keep.items():
            ui.update_selectize(
                id_, choices={"": f"Choose an enzyme… ({len(cutters)} cut once)", **{c: c for c in cutters}},
                selected=current if current in cutters else "",
            )

    def _flank(value: str):
        return value.strip() or None          # blank = work it out from the vector

    @reactive.calc
    def cloning_site():
        if input.enzyme_mode() == "two":
            pair = _two_enzymes()
            if pair is None:
                raise SilentException()
            return find_cloning_site(
                vector_record(), pair[0], input.orientation(), enzyme_b=pair[1],
                left_extra=_flank(input.left_extra()), right_extra=_flank(input.right_extra()),
            )
        return find_cloning_site(vector_record(), enzyme_name(), input.orientation())

    @reactive.effect
    def _show_suggested_flanks():
        """Show the flanking bases worked out from the vector as the placeholders."""
        pair = _two_enzymes() if input.enzyme_mode() == "two" else None
        if pair is None:
            return
        try:
            site = find_cloning_site(vector_record(), pair[0], input.orientation(), enzyme_b=pair[1])
        except Exception:
            return
        ui.update_text("left_extra", placeholder=f"auto: {site.left_extra or '(none)'}")
        ui.update_text("right_extra", placeholder=f"auto: {site.right_extra or '(none)'}")

    @reactive.calc
    def suggestions():
        """Cloning sites just downstream of the vector's Pol III promoter (may be empty)."""
        try:
            return suggest_cloning_sites(vector_record())
        except Exception:
            return []

    def _suggestion_hint(site=None):
        """A hint listing the sites that sit after the promoter, shown when the current one does not."""
        found = suggestions()
        if not found or (site is not None and site.promoter):
            return None
        names = "; ".join(x.enzymes for x in found[:6])
        return ui.p(_fa("lightbulb"), " Cloning sites downstream of the Pol III promoter in this vector: ",
                    ui.tags.strong(names), ". Choose one in step 2.", class_="hint small mt-2 mb-0")

    @render.ui
    def site_info():
        try:
            site = cloning_site()
        except SilentException:
            return ui.p((_missing_setup() or "Choose a vector and enzymes to see the cloning site."),
                        class_="hint mb-0")
        except Exception as e:
            return ui.TagList(_error_alert(e), _suggestion_hint())

        chips = [
            _chip("Enzyme" + ("s" if site.enzyme_b else ""), site.enzymes),
            _chip("Left overhang", f"{site.left_overhang} ({site.left_kind})"),
            _chip("Right overhang", f"{site.right_overhang} ({site.right_kind})"),
            _chip("Replaced", f"{site.removed_length} bp"),
            _chip("Orientation", site.orientation),
            _chip("Promoter", site.promoter or "none found upstream"),
        ]
        if site.left_extra or site.right_extra:
            chips += [_chip("Bases before spacer", site.left_extra or "–"),
                      _chip("Bases after spacer", site.right_extra or "–")]

        top, bottom = design_template(site, input.g_mode())

        def tpl_row(label, template):
            pre, spacer, post = split_spacer(template)
            return ui.div(
                ui.span(label, class_="tpl-label"),
                ui.tags.code("5′-", pre, ui.span(spacer, class_="tpl-spacer"), post, "-3′"),
                class_="tpl-row",
            )

        note = (" The G is added only when the sgRNA does not already start with one."
                if input.g_mode() == "if_needed" else "")
        return ui.TagList(
            ui.div(
                ui.div(_fa("circle-check", fill="var(--bs-success)"),
                       "The vector was digested in silico.", class_="mb-2 fw-semibold"),
                ui.div(*chips),
                ui.div(
                    tpl_row("Sense", top), tpl_row("Antisense", bottom),
                    ui.div(["Oligo design for a 20-nt spacer (", ui.span("N", class_="tpl-spacer"),
                            f" = your spacer; everything else comes from the vector).{note}"],
                           class_="hint small mt-1"),
                    class_="tpl",
                ),
                class_="status-ok ps-3",
            ),
            _suggestion_hint(site),
            *[_alert("warning", _fa("triangle-exclamation"), w) for w in site.warnings],
        )

    # --- Reactive calculations (run on Run / example clicks) ---

    @reactive.calc
    @reactive.event(input.run_btn, input.example_btn)
    def processed_rows() -> list[dict]:
        if _use_example():
            rows = read_input(_EXAMPLE_XLSX)
        else:
            req(input.sgrna_file())
            rows = read_input(Path(input.sgrna_file()[0]["datapath"]))
        validate_ids(rows, id_col=input.id_col())
        return add_cloning_oligos(
            rows, cloning_site(), seq_col=input.seq_col(), id_col=input.id_col(), g_mode=input.g_mode(),
        )

    @reactive.calc
    @reactive.event(input.run_btn, input.example_btn)
    def idt_rows() -> list[dict]:
        return build_idt_rows(processed_rows(), id_col=input.id_col())

    @reactive.calc
    @reactive.event(input.run_btn, input.example_btn)
    def plasmid_records():
        return build_plasmid_records(
            processed_rows(), vector_record(), cloning_site(),
            seq_col=input.seq_col(), id_col=input.id_col(),
        )

    def _all_warnings() -> list[str]:
        """Cloning-site warnings (e.g. unconfirmed orientation) plus every sgRNA's warnings."""
        return list(cloning_site().warnings) + [w for r in processed_rows() for w in r["warnings"]]

    # --- Notification on completion ---

    @reactive.effect
    @reactive.event(input.run_btn, input.example_btn)
    def _notify():
        try:
            n     = len(processed_rows())
            warns = _all_warnings()
            ui.notification_show(
                f"Designed {n * 2} oligos for {n} sgRNAs" + (f" ({len(warns)} warning(s))" if warns else ""),
                type="warning" if warns else "message", duration=4,
            )
        except SilentException:
            ui.notification_show(_missing_input() or "Something is missing from the inputs.",
                                 type="warning", duration=6)
        except Exception as e:
            ui.notification_show(str(e) if isinstance(e, ValueError) else f"Unexpected error: {e}",
                                 type="error", duration=8)

    # --- Results ---

    @render.ui
    def results_ui():
        if _stale():
            return None
        try:
            rows = processed_rows()
            warns = _all_warnings()
        except SilentException:
            if not _ran():
                return _empty_state()
            return _alert("warning", _fa("triangle-exclamation"),
                          _missing_input() or "Something is missing from the inputs.")
        except Exception as e:
            return _error_alert(e)

        n = len(rows)
        if n == 0:
            return _alert("warning", _fa("triangle-exclamation"), "The table contains no sgRNAs.")
        boxes = ui.layout_columns(
            ui.value_box("sgRNAs", f"{n}", showcase=icon_svg("dna", margin_left="0", margin_right="0"),
                         theme="primary", height="105px"),
            ui.value_box("Oligos to order", f"{2 * n}",
                         showcase=icon_svg("vial", margin_left="0", margin_right="0"), theme="primary",
                         height="105px"),
            ui.value_box("Warnings", f"{len(warns)}",
                         showcase=icon_svg("triangle-exclamation" if warns else "circle-check",
                                           margin_left="0", margin_right="0"),
                         theme="warning" if warns else "success", height="105px"),
            col_widths={"xs": 12, "md": 4}, fill=False,
        )
        warn_alert = _alert(
            "warning", ui.tags.strong(f"{len(warns)} warning(s)"),
            ui.tags.ul(*[ui.tags.li(w) for w in warns], class_="mb-0"),
        ) if warns else None
        downloads = ui.div(
            ui.download_button("download_all", "Everything (.zip)", icon=_fa("file-zipper"),
                               class_="btn-primary btn-sm"),
            ui.download_button("download_idt", "IDT order form", icon=_fa("file-excel"),
                               class_="btn-outline-primary btn-sm"),
            ui.download_button("download_gb", "GenBank (.gb)", icon=_fa("dna"),
                               class_="btn-outline-primary btn-sm"),
            ui.download_button("download_zip", "Per-sgRNA GenBank (.zip)", icon=_fa("file-zipper"),
                               class_="btn-outline-primary btn-sm"),
            ui.download_button("download_csv", "Oligos (.csv)", icon=_fa("file-csv"),
                               class_="btn-outline-primary btn-sm"),
            class_="header-actions",
        )
        return ui.TagList(
            boxes, warn_alert,
            ui.card(
                ui.card_header(_fa("table"), "Oligos", downloads),
                ui.output_data_frame("results_grid"),
                full_screen=True,
            ),
        )

    @render.data_frame
    def results_grid():
        req(not _stale())
        rows = processed_rows()
        df = pd.DataFrame(oligo_table(rows, id_col=input.id_col(), seq_col=input.seq_col()),
                          columns=TABLE_COLUMNS)
        df = df.drop(columns=["Spacer (5'→3')"]).rename(columns={"sgRNA": "ID"})
        flagged = [i for i, r in enumerate(rows) if r["warnings"]]
        styles = [{"cols": [1, 2], "style": {"font-family": "ui-monospace, Menlo, monospace", "font-size": "0.84rem"}}]
        if flagged:
            styles.append({"rows": flagged, "class": "table-warning"})
        return render.DataGrid(df, width="100%", summary=False, styles=styles)

    # --- Downloads with dynamic filenames ---

    def _base_name() -> str:
        if _use_example():
            return "example_sgRNA_seqs"
        f = input.sgrna_file()
        return sanitize_name(Path(f[0]["name"]).stem, "sgRNA") if f else "sgRNA"

    @render.download(filename=lambda: f"{_base_name()}_IDT_orderForm.xlsx")
    def download_idt():
        req(not _stale())
        req(idt_rows())
        with tempfile.TemporaryDirectory() as tmp_dir:
            path = Path(tmp_dir) / "order.xlsx"
            write_idt_excel(idt_rows(), path)
            yield path.read_bytes()

    @render.download(filename=lambda: f"{sanitize_name(input.combined_name(), 'combined_plasmids')}.gb")
    def download_gb():
        req(not _stale())
        req(plasmid_records())
        yield genbank_string([named_record(r, vector_label()) for r in plasmid_records()]).encode()

    @render.download(filename=lambda: f"{vector_label()}_plasmids.zip")
    def download_zip():
        req(not _stale())
        req(plasmid_records())
        yield genbank_zip(plasmid_records(), vector_label())

    @render.download(filename=lambda: f"{_base_name()}_oligos.csv")
    def download_csv():
        req(not _stale())
        req(processed_rows())
        yield oligos_csv(processed_rows(), id_col=input.id_col(), seq_col=input.seq_col())

    @render.download(filename=lambda: f"{_base_name()}_{vector_label()}_design.zip")
    def download_all():
        req(not _stale())
        req(plasmid_records())
        with tempfile.TemporaryDirectory() as tmp_dir:
            yield bundle_zip(
                idt_rows(), processed_rows(), plasmid_records(), vector_label(),
                sanitize_name(input.combined_name(), "combined_plasmids"), Path(tmp_dir) / "order.xlsx",
                id_col=input.id_col(), seq_col=input.seq_col(),
            )


# ---------------------------------------------------------------------------
# App
# ---------------------------------------------------------------------------

app = App(app_ui, server)

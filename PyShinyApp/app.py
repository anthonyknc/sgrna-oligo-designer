"""
app.py  —  sgRNA Oligo Designer v2.0.0
---------------------------------------
Core shiny mode. Run: shiny run app.py --reload
Deploy:        rsconnect deploy shiny . --name <account> --title sgrna-oligo-designer
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pandas as pd
import shinyswatch
from shiny import App, reactive, render, req, ui
from shiny.types import SilentException
from shiny import ui as sui

from sgrna_oligo_designer import (
    ENZYMES,
    add_cloning_oligos,
    all_enzymes,
    build_idt_rows,
    build_plasmid_records,
    find_cloning_site,
    read_input,
    read_vector_gb,
    sanitize_name,
    single_cutters,
    validate_ids,
    write_gb_files,
    write_idt_excel,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_VERSION      = "v2.0.0"
_LAST_UPDATE  = "2026-11-23"   # YYYY-MM-DD; bump together with _VERSION on each release
_GITHUB_URL   = "https://github.com/anthonyknc/sgrna-oligo-designer"
_HERE         = Path(__file__).parent
_EXAMPLE_XLSX = _HERE / "example_inputs" / "example_sgRNA_seqs.xlsx"

# Vectors bundled with the app: key -> (label, GenBank file)
_VECTORS = {
    "ipUSEPR":          ("ipUSEPR",          _HERE / "vectors" / "ipUSEPR.gb"),
    "lentiCRISPRv2":    ("lentiCRISPR v2",   _HERE / "vectors" / "lentiCRISPRv2.gb"),
    "LentiCRISPRv2GFP": ("LentiCRISPRv2GFP", _HERE / "vectors" / "LentiCRISPRv2GFP.gb"),
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
# Assets
# ---------------------------------------------------------------------------

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

# SVG icons for the dark mode button
_MOON_SVG = """
<svg xmlns="http://www.w3.org/2000/svg" width="16" height="16"
     viewBox="0 0 24 24" fill="currentColor">
  <path d="M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z"/>
</svg>"""

_SUN_SVG = """
<svg xmlns="http://www.w3.org/2000/svg" width="16" height="16"
     viewBox="0 0 24 24" fill="none" stroke="currentColor"
     stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
  <circle cx="12" cy="12" r="5"/>
  <line x1="12" y1="1"  x2="12" y2="3"/>
  <line x1="12" y1="21" x2="12" y2="23"/>
  <line x1="4.22"  y1="4.22"  x2="5.64"  y2="5.64"/>
  <line x1="18.36" y1="18.36" x2="19.78" y2="19.78"/>
  <line x1="1"  y1="12" x2="3"  y2="12"/>
  <line x1="21" y1="12" x2="23" y2="12"/>
  <line x1="4.22"  y1="19.78" x2="5.64"  y2="18.36"/>
  <line x1="18.36" y1="5.64"  x2="19.78" y2="4.22"/>
</svg>"""

_DARK_MODE_JS = """
Shiny.addCustomMessageHandler('set-theme', function(theme) {
    document.documentElement.setAttribute('data-bs-theme', theme);
    var moon = document.getElementById('dm-icon-moon');
    var sun  = document.getElementById('dm-icon-sun');
    if (moon && sun) {
        moon.style.display = (theme === 'dark') ? 'none'   : 'flex';
        sun.style.display  = (theme === 'dark') ? 'flex'   : 'none';
    }
});
"""

# CSS: spinner on Run button when Shiny is busy
_BUSY_CSS = """
#run_btn { transition: opacity 0.2s; position: relative; }
.shiny-busy #run_btn { opacity: 0.65; cursor: wait; }
.shiny-busy #run_btn::after {
    content: "";
    position: absolute;
    width: 0.9rem;
    height: 0.9rem;
    border: 2px solid rgba(255,255,255,0.7);
    border-top-color: transparent;
    border-radius: 50%;
    animation: shiny-spin 0.6s linear infinite;
    right: 0.7rem;
    top: 50%;
    transform: translateY(-50%);
}
@keyframes shiny-spin { to { transform: translateY(-50%) rotate(360deg); } }
"""

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _info(title: str, *content, placement: str = "right"):
    """
    A small '?' button that opens a popover when clicked (click only, so it also works on
    touch screens; click again or click elsewhere to close). A real <button>, so it is
    reachable by keyboard and does not trigger the surrounding label's input.
    """
    return sui.popover(
        sui.tags.button("?", type="button", class_="info-btn",
                        **{"aria-label": f"More information: {title}"}),
        *content,
        title=title,
        placement=placement,
        options={"trigger": "click"},
    )


def _label(text: str, *tip, title: str | None = None):
    """Input label followed by a click-to-open '?' popover."""
    return sui.tags.span(text, _info(title or text, *tip))


def _link(text: str, href: str):
    return sui.tags.a(text, href=href, target="_blank", rel="noopener")


_CRISPICK = "https://portals.broadinstitute.org/gppx/crispick/public"

_TIP_TABLE = (
    "A spreadsheet (.xlsx, .csv or .tsv) with one row per sgRNA and two columns: ",
    sui.tags.code("sgRNA_id"), " (unique, no spaces) and ", sui.tags.code("sgRNA_seq"),
    " (the 20-nt target sequence, without the PAM or any overhangs). Extra columns are ignored, "
    "and spaces and lowercase are cleaned automatically. Design guides with ",
    _link("CRISPick", _CRISPICK), ". More under Guide → Getting started.",
)
_TIP_VECTOR = (
    "The empty vector you will clone into. ipUSEPR, lentiCRISPR v2 and LentiCRISPRv2GFP are "
    "built in. Choose “Upload my own” for another circular GenBank file.",
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
    return sui.markdown((_HERE / "help" / f"{name}.md").read_text(encoding="utf-8"))


# Close any open info popover on an outside click or Escape. (Bootstrap's click trigger only
# toggles; its 'focus' trigger would break in Safari, which doesn't focus buttons on click.)
_INFO_JS = """
(function () {
  function closeOpen(except) {
    document.querySelectorAll('.popover.show').forEach(function (pop) {
      var trigger = document.querySelector('[aria-describedby="' + pop.id + '"]');
      if (!trigger || trigger === except) return;
      var inst = window.bootstrap && bootstrap.Popover.getInstance(trigger);
      if (inst) inst.hide();
    });
  }
  document.addEventListener('click', function (e) {
    if (e.target.closest('.popover')) return;
    closeOpen(e.target.closest('.info-btn'));
  });
  document.addEventListener('keydown', function (e) { if (e.key === 'Escape') closeOpen(null); });
})();
"""

_INFO_CSS = """
.info-btn {
    display: inline-flex; align-items: center; justify-content: center;
    width: 1.15rem; height: 1.15rem; margin-left: 0.4rem; padding: 0;
    border: 1px solid #adb5bd; border-radius: 50%; background: transparent;
    color: #6c757d; font-size: 0.7rem; font-weight: 700; line-height: 1;
    cursor: pointer; vertical-align: middle;
}
.info-btn:hover, .info-btn:focus-visible { color: #fff; background: #6c757d; border-color: #6c757d; }
.popover { max-width: 320px; }
.guide h1 { font-size: 1.6rem; margin-bottom: 1rem; }
.guide h2 { font-size: 1.25rem; margin-top: 1.75rem; }
.guide { max-width: 52rem; }
"""


# ---------------------------------------------------------------------------
# Title bar
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

        # --- Load example ---
        sui.input_action_button(
            "example_btn",
            "⚡ Load example sgRNAs",
            class_="btn-outline-secondary w-100",
        ),

        sui.hr(),

        # --- Step 1 ---
        sui.h6("Step 1 — sgRNA table"),
        sui.input_file(
            "sgrna_file",
            _label("Upload sgRNA table", *_TIP_TABLE, title="sgRNA table"),
            accept=[".xlsx", ".xlsm", ".csv", ".tsv"],
        ),

        sui.hr(),

        # --- Step 2 ---
        sui.h6("Step 2 — Vector and enzyme"),
        sui.input_select(
            "vector_choice",
            _label("Vector", *_TIP_VECTOR),
            {**{k: v[0] for k, v in _VECTORS.items()}, _UPLOAD: "Upload my own vector (.gb)…"},
        ),
        sui.panel_conditional(
            f"input.vector_choice == '{_UPLOAD}'",
            sui.input_file(
                "vector_file",
                "Upload empty vector (.gb / .gbk)",
                accept=[".gb", ".gbk"],
            ),
        ),
        sui.input_radio_buttons(
            "enzyme_mode",
            _label("Cloning enzymes", *_TIP_MODE),
            {"one": "One enzyme (cuts the vector twice)", "two": "Two different enzymes (each cuts once)"},
            selected="one",
        ),
        sui.panel_conditional(
            "input.enzyme_mode == 'one'",
            sui.input_select(
                "enzyme",
                _label("Restriction enzyme", *_TIP_ENZYME),
                {**{e: _ENZYME_LABELS.get(e, e) for e in ENZYMES}, _OTHER: "Other…"},
            ),
            sui.panel_conditional(
                f"input.enzyme == '{_OTHER}'",
                sui.input_selectize(
                    "other_enzyme",
                    "Search all enzymes",
                    {n: f"{n}  ({e.site}, {e.overhang_len}-nt overhang)" for n, e in all_enzymes().items()},
                    selected="SapI",
                ),
            ),
        ),
        sui.panel_conditional(
            "input.enzyme_mode == 'two'",
            sui.input_selectize("enzyme_a", "Enzyme 1 (cuts once)", {"": "Choose an enzyme…"}),
            sui.input_selectize("enzyme_b", "Enzyme 2 (cuts once)", {"": "Choose an enzyme…"}),
        ),
        sui.input_text(
            "vector_name",
            _label("Vector name (optional)",
                   "Name prepended to each plasmid name in the GenBank file. "
                   "Defaults to the vector name or the uploaded file's name."),
            placeholder="Default: vector name",
        ),
        sui.input_text(
            "combined_name",
            _label("Combined .gb filename stem",
                   "The combined GenBank file will be named <stem>.gb."),
            value="combined_plasmids",
        ),

        sui.hr(),

        # --- Advanced options (collapsible) ---
        sui.accordion(
            sui.accordion_panel(
                "⚙ Advanced options",
                sui.input_text(
                    "seq_col",
                    _label("sgRNA sequence column",
                           "Column name in your input file containing the 20-nt sgRNA target sequences."),
                    value="sgRNA_seq",
                ),
                sui.input_text(
                    "id_col",
                    _label("sgRNA ID column",
                           "Column name in your input file containing unique sgRNA identifiers."),
                    value="sgRNA_id",
                ),
                sui.panel_conditional(
                    "input.enzyme_mode == 'two'",
                    sui.input_text("left_extra", _label("Bases before the spacer", *_TIP_FLANK, title="Flanking bases"),
                                   placeholder="auto"),
                    sui.input_text("right_extra", _label("Bases after the spacer", *_TIP_FLANK, title="Flanking bases"),
                                   placeholder="auto"),
                ),
                sui.input_radio_buttons(
                    "g_mode",
                    _label("5' G for U6 transcription", *_TIP_G),
                    _G_CHOICES,
                    selected="always",
                ),
                sui.input_select(
                    "orientation",
                    _label("sgRNA orientation in the vector",
                           "Which strand of the GenBank file the sgRNA is transcribed from. "
                           "Auto-detect looks for the start of the Cas9 scaffold (GTTT)."),
                    _ORIENT_CHOICES,
                ),
            ),
            open=False,
        ),

        sui.hr(),
        sui.input_action_button("run_btn", "Run", class_="btn-primary w-100"),

        width=320,
    ),

    # --- Dark mode toggle: circular icon button, top-right of navbar ---
    sui.tags.div(
        sui.input_action_button(
            "dark_mode_btn",
            sui.HTML(
                f'<span id="dm-icon-moon" style="display:flex;">'
                f'{_MOON_SVG}</span>'
                f'<span id="dm-icon-sun" style="display:none;">'
                f'{_SUN_SVG}</span>'
            ),
            title="Toggle dark / light mode",
            style=(
                "width:2.2rem; height:2.2rem; padding:0; border-radius:50%;"
                " background:transparent;"
                " border:1px solid rgba(128,128,128,0.3);"
                " cursor:pointer; display:inline-flex;"
                " align-items:center; justify-content:center;"
                " color:inherit; transition:border-color 0.2s, opacity 0.2s;"
            ),
        ),
        style="position:fixed; top:0.5rem; right:0.75rem; z-index:1031;",
    ),

    # --- JS & CSS ---
    sui.tags.script(_DARK_MODE_JS),
    sui.tags.script(_INFO_JS),
    sui.tags.style(_BUSY_CSS + _INFO_CSS),

    # --- Main panel: Designer + Guide tabs ---
    sui.navset_tab(
        sui.nav_panel("Designer", sui.tags.div(
        sui.h4("Cloning site", _info("Cloning site", *_TIP_SITE, placement="bottom")),
        sui.output_ui("site_info"),
        sui.output_ui("stale_notice"),

        sui.h4("Summary"),
        sui.output_text("summary"),
        sui.output_ui("warnings_ui"),

        sui.hr(),
        sui.h4("IDT order form preview"),
        sui.output_table("idt_table"),

        sui.hr(),
        sui.h4("Downloads"),
        sui.h6("IDT Bulk Order Form (.xlsx)"),
        sui.download_button("download_idt", "Download IDT Order Form"),
        sui.h6("Combined GenBank File (.gb)"),
        sui.download_button("download_gb", "Download Combined GenBank")
        , style="padding-top:1rem;")),
        sui.nav_panel(
            "Guide",
            sui.tags.div(
                sui.navset_pill(
                    sui.nav_panel("Getting started", sui.tags.div(_guide("getting_started"), class_="guide")),
                    sui.nav_panel("How it works", sui.tags.div(_guide("how_it_works"), class_="guide")),
                    sui.nav_panel("Troubleshooting", sui.tags.div(_guide("troubleshooting"), class_="guide")),
                    sui.nav_panel("Citing", sui.tags.div(_guide("citing"), class_="guide")),
                ),
                style="padding-top:1rem;",
            ),
        ),
    ),

    # --- Footer ---
    sui.tags.div(
        sui.tags.hr(),
        sui.tags.p(
            "© 2026 Anthony Chan · ",
            sui.tags.a("GitHub", href=_GITHUB_URL, target="_blank"),
            " · MIT License",
            style="text-align:center; color:#6c757d; font-size:0.8em; margin:0;",
        ),
        sui.tags.p(
            f"Last Update: {_LAST_UPDATE}",
            style="text-align:center; color:#6c757d; font-size:0.8em; margin:0.25rem 0 0;",
        ),
        style="margin-top:2rem; padding-bottom:1rem;",
    ),

    title=_TITLE,
    theme=shinyswatch.theme.flatly(),
)


# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

def server(input, output, session):

    # Example mode flag
    _use_example = reactive.value(False)

    @reactive.effect
    @reactive.event(input.example_btn)
    def _set_example():
        _use_example.set(True)

    @reactive.effect
    @reactive.event(input.run_btn)
    def _clear_example():
        _use_example.set(False)

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

    @output
    @render.ui
    def stale_notice():
        if _stale():
            return sui.tags.div(
                "Settings changed since the last run. Click Run to update the results.",
                class_="alert alert-info", role="alert",
            )

    # --- Dark mode ---

    _dark = reactive.value(False)

    @reactive.effect
    @reactive.event(input.dark_mode_btn)
    async def _apply_theme():
        if input.dark_mode_btn() > 0:   # ignore initial value of 0
            new_state = not _dark()
            _dark.set(new_state)
            await session.send_custom_message(
                "set-theme", "dark" if new_state else "light"
            )

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

    @reactive.effect
    def _offer_single_cutters():
        """In two-enzyme mode, offer only the enzymes that cut this vector exactly once."""
        if input.enzyme_mode() != "two":
            return
        try:
            cutters = single_cutters(vector_record())
        except Exception:
            return
        with reactive.isolate():
            keep = {"enzyme_a": input.enzyme_a(), "enzyme_b": input.enzyme_b()}
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

    @output
    @render.ui
    def site_info():
        try:
            site = cloning_site()
        except SilentException:
            return sui.tags.p("Choose the two enzymes to see the cloning site.", style="color:#6c757d;")
        except Exception as e:
            if input.vector_choice() == _UPLOAD and not input.vector_file():
                return sui.tags.p("Upload a vector GenBank file to see its cloning site.",
                                  style="color:#6c757d;")
            return sui.tags.div(str(e), class_="alert alert-danger", role="alert")
        cuts = "cuts the vector once each" if site.enzyme_b else "cuts the vector at two sites"
        flanks = []
        if site.left_extra or site.right_extra:
            flanks = [
                " Constant bases added around the spacer: ",
                sui.tags.code(site.left_extra or "–"), " before and ",
                sui.tags.code(site.right_extra or "–"), " after.",
            ]
        items = [
            sui.tags.div(
                sui.tags.strong(f"✓ {site.enzymes}"), f" {cuts}. ",
                "Overhangs: ", sui.tags.code(site.left_overhang), f" ({site.left_kind} overhang, left) and ",
                sui.tags.code(site.right_overhang), f" ({site.right_kind} overhang, right). ",
                f"{site.removed_length} bp between the cuts are replaced by the sgRNA. "
                f"Orientation: {site.orientation}.", *flanks,
                class_="alert alert-success", role="alert",
            )
        ]
        items += [sui.tags.div(w, class_="alert alert-warning", role="alert") for w in site.warnings]
        return sui.TagList(*items)

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
            rows,
            cloning_site(),
            seq_col=input.seq_col(),
            id_col=input.id_col(),
            g_mode=input.g_mode(),
        )

    @reactive.calc
    @reactive.event(input.run_btn, input.example_btn)
    def idt_rows() -> list[dict]:
        return build_idt_rows(processed_rows(), id_col=input.id_col())

    @reactive.calc
    @reactive.event(input.run_btn, input.example_btn)
    def plasmid_records():
        return build_plasmid_records(
            processed_rows(),
            vector_record(),
            cloning_site(),
            seq_col=input.seq_col(),
            id_col=input.id_col(),
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
                f"✓ Processed {n} sgRNA sequences → {n * 2} oligos"
                + (f" ({len(warns)} warning(s), see below)" if warns else ""),
                type="warning" if warns else "message",
                duration=4,
            )
        except SilentException:
            ui.notification_show(
                "Upload an sgRNA table (or click Load example sgRNAs), and a vector file if you chose "
                "'Upload my own', then click Run.",
                type="warning", duration=6,
            )
        except Exception as e:
            ui.notification_show(str(e), type="error", duration=8)

    # --- Outputs ---

    @output
    @render.text
    def summary():
        req(not _stale())
        req(processed_rows())
        n = len(processed_rows())
        return (
            f"Processed {n} sgRNA sequences "
            f"→ {n * 2} oligos ({n} sense + {n} antisense) | GenBank file ready"
        )

    @output
    @render.ui
    def warnings_ui():
        req(not _stale())
        req(processed_rows())
        warns = _all_warnings()
        if not warns:
            return None
        return sui.tags.div(
            sui.tags.strong(f"{len(warns)} warning(s)"),
            sui.tags.ul(*[sui.tags.li(w) for w in warns]),
            class_="alert alert-warning", role="alert",
        )

    @output
    @render.table
    def idt_table():
        req(not _stale())
        req(idt_rows())
        return pd.DataFrame(idt_rows())

    # --- Downloads with dynamic filenames ---

    @render.download(
        filename=lambda: (
            "example_sgRNA_seqs_IDT_orderForm.xlsx"
            if _use_example()
            else f"{sanitize_name(Path(input.sgrna_file()[0]['name']).stem, 'sgRNA')}_IDT_orderForm.xlsx"
            if input.sgrna_file()
            else "IDT_orderForm.xlsx"
        )
    )
    def download_idt():
        req(not _stale())
        req(idt_rows())
        with tempfile.NamedTemporaryFile(suffix=".xlsx", delete=False) as tmp:
            tmp_path = Path(tmp.name)
        try:
            write_idt_excel(idt_rows(), tmp_path)
            yield tmp_path.read_bytes()
        finally:
            tmp_path.unlink(missing_ok=True)

    @render.download(
        filename=lambda: f"{sanitize_name(input.combined_name(), 'combined_plasmids')}.gb"
    )
    def download_gb():
        req(not _stale())
        req(plasmid_records())
        vector_name  = vector_label()
        combined_stem = sanitize_name(input.combined_name(), "combined_plasmids")
        with tempfile.TemporaryDirectory() as tmp_dir:
            write_gb_files(
                plasmid_records(),
                output_dir=Path(tmp_dir),
                vector_name=vector_name,
                multi_gb=True,
                combined_stem=combined_stem,
            )
            yield (Path(tmp_dir) / f"{combined_stem}.gb").read_bytes()


# ---------------------------------------------------------------------------
# App
# ---------------------------------------------------------------------------

app = App(app_ui, server)

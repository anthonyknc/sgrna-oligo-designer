# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

# sgRNA Oligo Designer

Designs sgRNA cloning oligos (IDT bulk order `.xlsx`) and annotated plasmid maps (`.gb`) for
Type IIS (restriction-enzyme) cloning of SpCas9 sgRNAs into lentiviral vectors. Audience: wet-lab scientists.

## Layout

- `PyShinyApp/` is the self-contained deploy dir for shinyapps.io: `app.py` (Shiny for Python,
  Core mode), `sgrna_oligo_designer.py` (core logic, no CLI), `help/*.md` (SOURCE OF TRUTH for the
  Guide tab; the README links to these files, keep their links absolute), `vectors/` (the three bundled
  vectors, offered as presets), `example_inputs/` (table for the "Load example data" button),
  and `requirements.txt` (GENERATED, don't hand-edit).
- `example_inputs/sgRNA_seqs/` and `example_outputs/` at the root are the reference data the tests use.
- `tests/`: `test_golden.py` (example vectors vs `example_outputs/`), `test_cloning.py` (synthetic
  vectors: every enzyme, both layouts, both orientations, errors), `test_core.py` (validation, I/O).
  `pythonpath` is set to `PyShinyApp` in `pyproject.toml`.
- uv project: `pyproject.toml` + `uv.lock`, Python 3.12 (`.python-version`, pinned
  `>=3.12,<3.13` until a test deploy says otherwise).

## Releasing

On each release bump `_VERSION` and `_LAST_UPDATE` in `PyShinyApp/app.py` and `version` in
`pyproject.toml` (a test checks that the first and last agree), then `uv lock`.

## Commands

```bash
uv sync                                        # create/refresh .venv from uv.lock
uv run pytest                                  # all tests
uv run pytest tests/test_core.py::test_default_oligos   # one test (add -k <expr> to filter)
uv run shiny run PyShinyApp/app.py --reload    # run the app locally
uv export --no-dev --no-hashes --no-emit-project -o PyShinyApp/requirements.txt   # after changing deps
uv run rsconnect deploy shiny PyShinyApp --name <account> --title sgrna-oligo-designer
```

No linter or formatter is configured.

## Architecture

`app.py` is a thin reactive layer over `sgrna_oligo_designer.py`; all logic lives in the core.
Rows are plain `list[dict]`, mutated in place.

Pipeline: `find_cloning_site(vector, enzyme, orientation)` digests the vector in silico and returns
a frozen `CloningSite` (overhangs, the region to replace, orientation, junction context). Everything
downstream is derived from it: `add_cloning_oligos(rows, site)` -> `build_idt_rows` ->
`write_idt_excel`, and `build_plasmid_records(rows, vector, site)` -> `insert_sgrna` ->
`write_gb_files`.

Things that span functions and are easy to miss:
- All `CloningSite` coordinates are in the *sense frame*. For a vector whose sgRNA is transcribed
  from the reverse strand, `insert_sgrna` reverse-complements the record, edits, and reverse-
  complements back; never mix frames.
- Oligos are `left_overhang + core` and `rc(core + right_overhang)` where `core` is the spacer with
  the optional 5' G (`g_mode`). The extra G/C of the classic `CACCG`/`AAAC...C` design fall out of
  this; they are not hard-coded.
- The stuffer is "the shorter circular piece between the two cuts"; both layouts (sites lost in the
  stuffer vs left in the backbone) use the same arithmetic. Cut geometry comes from `Bio.Restriction`.
- If the vector already has the 5' G after the left overhang (ipUSEPR, lentiCRISPR), `insert_sgrna`
  keeps that base and its annotations and replaces only what follows.
- Warnings are returned in `row["warnings"]` / `CloningSite.warnings` (stderr is invisible in the
  web app). The app shows them.
- Field help is a click-only `?` popover (`_info`/`_label` in `app.py`); a small JS handler closes
  open popovers on outside click/Escape. Don't switch to Bootstrap's `focus` trigger: Safari doesn't
  focus buttons on click.
- User-supplied text that reaches files is sanitised: sgRNA IDs go through `validate_ids` (letters,
  digits, `_ - .` only; blocks spaces that break GenBank LOCUS lines and a leading `=` that Excel
  runs as a formula) and file/vector names through `sanitize_name` (cannot leave the output dir).
  `write_idt_excel` also forces string cells. Keep any new user-controlled name on this path.
- The app's `cloning_site` calc is live (updates on vector/enzyme change and shows errors inline);
  results calcs fire only on Run/example clicks and are hidden as stale when a setting changes.
- Golden tests compare GenBank output byte-for-byte (including the LOCUS line), so any change to
  record naming, feature order or qualifiers shows up as a golden diff. The `source` feature in the
  goldens spans the whole new plasmid (an intentional correction of the old stuffer-based output).

## Decisions (made with the user; don't re-litigate)

- **Shiny app is the only entry point.** The CLI is removed. Keep the core as a plain
  importable module with tests; the app stays thin.
- **pixi -> uv.** Single `pyproject.toml` + `uv.lock`; `requirements.txt` for shinyapps.io is
  generated via `uv export`. `rsconnect-python` is a dev dependency, not a runtime one.
- **Stay on Biopython.** biotite was considered and rejected (no real speed benefit for ~10 kb
  files, GenBank round-trip fidelity with SnapGene/Benchling matters, and `Bio.Restriction`
  is needed for the enzyme work).
- **Keep the current Python version (3.12, see `PyShinyApp/.python-version`)** for now. Test a
  shinyapps.io deploy early before changing anything about the Python version.
- **Enzyme-based design (DONE; replaced the manual overhang flags and the stuffer input).** User picks vector + enzyme; the tool
  finds the flanking sites, simulates the cut, derives the 4-nt overhangs, and builds oligos and
  the plasmid from the real cut ends. The stuffer is inferred (fragment between cuts), not typed.
  - Enzymes: BsmBI, Esp3I, BsaI, BbsI (`ENZYMES`; add more by name, geometry comes from Bio.Restriction).
  - Support **both** layouts: recognition sites inside the stuffer (lost on cutting, e.g.
    lentiCRISPR) and sites that stay in the backbone.
  - Warn if a spacer contains the enzyme site on either strand; error if the vector has
    unexpected extra sites.
  - The 5' G for U6 is a single option, not derived from the vector.
- Test vectors: ipUSEPR first, then LentiCRISPRv2 and LentiCRISPRv2GFP from `example_inputs/`.

## Work order

1. DONE: golden tests.
2. DONE: uv migration + restructure.
3. Mostly done: `uv lock` resolved the latest compatible versions. Remaining: a real shinyapps.io test
   deploy (needs the user's credentials; do not run it unprompted).
4. DONE: enzyme-based redesign. Only BsmBI vectors exist in the repo; BsaI/BbsI/Esp3I are covered by
   synthetic-vector tests, not by a real vector. Worth testing with a real BsaI vector if one turns up.

## Known issues

- The app only offers the combined `.gb`; per-sgRNA files are not offered.
- A cloning site that spans the plasmid origin is rejected rather than handled.
- A local git `post-commit` hook used to copy the root script into `PyShinyApp/`; it is disabled
  (`.git/hooks/post-commit.disabled`) because there is now a single copy.

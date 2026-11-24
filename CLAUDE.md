# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

# sgRNA Oligo Designer

Designs sgRNA cloning oligos (IDT bulk order `.xlsx`) and annotated plasmid maps (`.gb`) for
Golden Gate cloning of SpCas9 sgRNAs into lentiviral vectors. Audience: wet-lab scientists.

## Layout

- `PyShinyApp/` is the self-contained deploy dir for shinyapps.io: `app.py` (Shiny for Python,
  Core mode), `sgrna_oligo_designer.py` (core logic, no CLI), `example_inputs/` (bundled for the
  "Load example data" button), and `requirements.txt` (GENERATED, don't hand-edit).
- `example_inputs/` and `example_outputs/` at the root are the reference data the tests use.
- `tests/`: golden tests against `example_outputs/` plus unit tests (`pythonpath` is set to
  `PyShinyApp` in `pyproject.toml`).
- uv project: `pyproject.toml` + `uv.lock`, Python 3.12 (`.python-version`, pinned
  `>=3.12,<3.13` until a test deploy says otherwise).

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

`app.py` is a thin reactive layer over `sgrna_oligo_designer.py`; all logic lives in the core and
the app only wires inputs to it. Pipeline (rows are plain `list[dict]`, mutated in place):

`read_input` -> `check_unique_ids` -> `add_cloning_oligos` (validates/cleans each sgRNA, adds
`oligo_sense`/`oligo_antisense`) -> `build_idt_rows` -> `write_idt_excel`; and, when a vector is
supplied, `read_vector_gb` -> `build_plasmid_records` -> `insert_sgrna` (uses `validate_stuffer` to
find the stuffer on either strand, splices in the sgRNA, drops features inside the stuffer, shifts
downstream features) -> `write_gb_files`.

Things that span files and are easy to miss:
- The sense oligo always gets `CACCG` (extra G for U6), but the plasmid insert is the bare 20-mer;
  the plasmid sequence is built by replacing the stuffer, not by simulating the digestion.
- `write_gb_files` renames records to `<vector_name>_<sgRNA_id>`; the app hard-codes `ipUSEPR` as
  the vector name in example mode.
- The app's reactive calcs fire only on Run/example clicks (`@reactive.event`), and clicking Run
  clears example mode.
- Golden tests compare GenBank output byte-for-byte (including the LOCUS line), so any change to
  record naming, feature order or qualifiers shows up as a golden diff.

## Decisions (made with the user; don't re-litigate)

- **Shiny app is the only entry point.** The CLI is being removed. Keep the core as a plain
  importable module with tests; the app stays thin.
- **pixi -> uv.** Single `pyproject.toml` + `uv.lock`; `requirements.txt` for shinyapps.io is
  generated via `uv export`. `rsconnect-python` is a dev dependency, not a runtime one.
- **Stay on Biopython.** biotite was considered and rejected (no real speed benefit for ~10 kb
  files, GenBank round-trip fidelity with SnapGene/Benchling matters, and `Bio.Restriction`
  is needed for the enzyme work).
- **Keep the current Python version (3.12, see `PyShinyApp/.python-version`)** for now. Test a
  shinyapps.io deploy early before changing anything about the Python version.
- **Enzyme-based design (replaces manual overhang flags).** User picks vector + enzyme; the tool
  finds the flanking sites, simulates the cut, derives the 4-nt overhangs, and builds oligos and
  the plasmid from the real cut ends. The stuffer is inferred (fragment between cuts), not typed.
  - Enzymes at launch: **BsmBI and BsaI** at minimum.
  - Support **both** layouts: recognition sites inside the stuffer (lost on cutting, e.g.
    lentiCRISPR) and sites that stay in the backbone.
  - Warn if a spacer contains the enzyme site on either strand; error if the vector has
    unexpected extra sites.
  - The 5' G for U6 is a single option, not derived from the vector.
- Test vectors: ipUSEPR first, then LentiCRISPRv2 and LentiCRISPRv2GFP from `example_inputs/`.

## Work order

1. DONE: golden tests against `example_outputs/`.
2. DONE: uv migration + restructure (CLI and pixi removed, single core module). Note `uv lock`
   resolved the latest compatible versions, so step 3 is mostly done; remaining: test deploy.
3. Upgrade all dependencies to latest compatible versions (smoke-test the app after).
4. Enzyme-based redesign.
5. Re-run `/init` after step 2 to refresh commands/architecture here.

## Known issues to fix during the work

- `insert_sgrna` only shifts features that start after the stuffer; features spanning the
  insertion site keep a stale end coordinate when sgRNA and stuffer lengths differ (they usually
  do: ipUSEPR stuffer is 17 nt, lentiCRISPRv2 is 1,880 bp). Covered by a strict `xfail` test.
- The app only offers the combined `.gb`; per-sgRNA files are not offered.
- A local git `post-commit` hook used to copy the root script into `PyShinyApp/`; it is disabled
  (`.git/hooks/post-commit.disabled`) because there is now a single copy.

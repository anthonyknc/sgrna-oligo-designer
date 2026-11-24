# sgRNA Oligo Designer

Designs sgRNA cloning oligos (IDT bulk order `.xlsx`) and annotated plasmid maps (`.gb`) for
Golden Gate cloning of SpCas9 sgRNAs into lentiviral vectors. Audience: wet-lab scientists.

## Current state (transitional)

- `sgrna_oligo_designer.py` (root) is the core + CLI. `PyShinyApp/sgrna_oligo_designer.py` is a
  byte-identical copy used by `PyShinyApp/app.py` (Shiny for Python, deployed to shinyapps.io).
  Until the restructure lands, any core change must be made in both copies.
- Environment is currently pixi (`pixi.toml`); moving to uv.

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

1. Golden tests against `example_outputs/` (in progress; guards the later steps).
2. uv migration + restructure into one core module + app.
3. Upgrade all dependencies to latest compatible versions (smoke-test the app after).
4. Enzyme-based redesign.
5. Re-run `/init` after step 2 to refresh commands/architecture here.

## Known issues to fix during the work

- `insert_sgrna` only shifts features that start after the stuffer; features spanning the
  insertion site keep a stale end coordinate when sgRNA and stuffer lengths differ (they usually
  do: ipUSEPR stuffer is 17 nt, lentiCRISPRv2 is 1,880 bp).
- App skips the duplicate-ID check that the CLI does, and doesn't pass `id_col` to
  `add_cloning_oligos`.
- Stale docs: module docstring flags (`--top-prefix`), IDT names are `_s`/`_a` not `_top`/`_bottom`,
  `--vector-name` help has the filename order backwards.

## Running tests

No dev tooling is in the repo yet. From the repo root:

```bash
uv run --python 3.12 --with pytest --with biopython --with openpyxl pytest
```

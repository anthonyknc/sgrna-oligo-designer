# Change log

Versions are `major.minor`. The version in the app header, `_VERSION` in `PyShinyApp/app.py` and
`version` in `pyproject.toml` move together (see CLAUDE.md, "Releasing"). Newest first.

## v2.1 (2026-11-24)

- Pol III promoter detection (human U6, mouse U6, H1; 7SK from annotations only): by annotation first,
  then by sequence. The cloning site shows a Promoter chip, and the promoter helps decide the orientation.
- Choosing a vector now selects its enzyme, or enzyme pair, automatically.
- Nine bundled vectors, each labelled with its Addgene number (or PMID) and default enzyme(s): new are
  lentiGuide-Puro, pX330, pU6-(BbsI)_CBh-Cas9-T2A-mCherry-H1, phH1-gRNA, pU6-sgRNA EF1Alpha-puro-T2A-BFP
  and pMCB306. lentiCRISPR v2 and LentiCRISPRv2GFP were replaced by the Addgene records (same sequences).
- New interface: navbar with Design and Guide tabs, light and dark mode, help popovers, including one for
  Load example sgRNAs.
- Header: GitHub and light/dark icons aligned; version shown as `vMAJOR.MINOR`.
- Two-enzyme cloning (for example BstXI + BlpI), with the flanking bases suggested from the scaffold.
- Enzyme-derived overhangs and stuffer; the Shiny app is the only entry point.

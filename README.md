# sgRNA Oligo Designer

A web app for designing sgRNA cloning oligos for CRISPR experiments.
Given a list of sgRNA sequences, it generates:

1. **An IDT bulk order form** (`.xlsx`) — ready to upload directly to [IDT](https://www.idtdna.com) for oligo synthesis
2. **Annotated plasmid maps** (`.gb`) — one GenBank record per sgRNA construct, in a single combined file, ready to open in [SnapGene](https://www.snapgene.com) or [Benchling](https://www.benchling.com) for sequencing verification

You choose a vector and a Type IIS restriction enzyme (*Bsm*BI/Esp3I, *Bbs*I/BpiI, *Bsa*I and SapI are listed first; any other commercially available Type IIS enzyme can be chosen under **Other…**, and vectors cut by two different enzymes, such as BstXI + BlpI, are supported too). The app cuts the
vector in silico and reads the overhangs and the stuffer from it, so there is nothing to type in by hand.

> 🌐 **https://anthonyknc.shinyapps.io/sgrna-oligo-designer**

No installation or coding required. Click **Load example sgRNAs** in the app to try it first. The app has a light and a dark mode, and the **Guide** tab explains every option.

---

## Guide

The same guide is built into the app (**Guide** tab, and the **?** icons next to each input).
These files are the single source of truth for it:

- [Getting started](PyShinyApp/help/getting_started.md): designing sgRNAs with CRISPick, the input table, running the app, all inputs
- [How it works](PyShinyApp/help/how_it_works.md): how the cloning is simulated, why there is an extra G, the outputs, supported vectors
- [Troubleshooting](PyShinyApp/help/troubleshooting.md): what each error and warning means
- [Citing](PyShinyApp/help/citing.md): how to cite this tool and the ipUSEPR vector

Example inputs are in [`example_inputs/`](example_inputs/), the bundled vectors are in
[`PyShinyApp/vectors/`](PyShinyApp/vectors/), and example outputs are in [`example_outputs/`](example_outputs/).

---

## Development

The app is a Shiny for Python application in [`PyShinyApp/`](PyShinyApp/), with the core
logic in `PyShinyApp/sgrna_oligo_designer.py`. Dependencies are managed with
[uv](https://docs.astral.sh/uv/) (Python 3.12).

```bash
git clone https://github.com/anthonyknc/sgrna-oligo-designer.git
cd sgrna-oligo-designer
uv sync                                       # create .venv from uv.lock
uv run shiny run PyShinyApp/app.py --reload   # run locally
uv run pytest                                 # run the tests
```

---

## Citing

If you use this tool in your research, please cite it as:

> Chan, A.K.N. (2026). sgRNA Oligo Designer. GitHub: https://github.com/anthonyknc/sgrna-oligo-designer.

More citations (the ipUSEPR vector, STAgR dual sgRNAs) are in [Citing](PyShinyApp/help/citing.md).

---

## License

MIT License. See [LICENSE](LICENSE) for details.

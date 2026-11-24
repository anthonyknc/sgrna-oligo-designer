# sgRNA Oligo Designer

A web app for designing sgRNA cloning oligos for CRISPR experiments.
Given a list of sgRNA sequences, it generates:

1. **An IDT bulk order form** (`.xlsx`) — ready to upload directly to [IDT](https://www.idtdna.com) for oligo synthesis
2. **Annotated plasmid maps** (`.gb`) — one GenBank record per sgRNA construct, in a single combined file, ready to open in [SnapGene](https://www.snapgene.com) or [Benchling](https://www.benchling.com) for sequencing verification

---

## Using the app

> 🌐 **https://anthonyknc.shinyapps.io/sgrna-oligo-designer**

No installation or coding required. Upload your sgRNA table, provide your vector GenBank file and stuffer sequence, then download your IDT order form and annotated GenBank file directly in the browser. Click **⚡ Load example data** to try it with the ipUSEPR example first.

---

## Step 1 — Design your sgRNA sequences with CRISPRpick

Before using this tool, you need a list of 20-nucleotide sgRNA target sequences for your gene(s) of interest. We recommend **CRISPRpick** from Broad Institute:

> 🔗 https://portals.broadinstitute.org/gppx/crispick/public

CRISPRpick is free to use and requires no account. Simply enter your gene name or transcript ID, select your genome and Cas9 system (e.g. *Sp*Cas9), and the tool will return ranked sgRNA sequences with on-target efficiency scores. Select your top candidates (typically 3–6 per gene) and export them to a spreadsheet, then reformat the exported data to match the input format shown in the examples at [`example_inputs/sgRNA_seqs/`](example_inputs/sgRNA_seqs/).

---

## Step 2 — Prepare your input file

Prepare a spreadsheet (Excel `.xlsx`, `.csv`, or `.tsv`) with at least these two columns:

| sgRNA_id | sgRNA_seq |
|---|---|
| sg_Luc_474 | GATTCTAAAACGGATTACCA |
| sg_tRFP657_294 | GTCACCACATACGAAGACGG |

- **`sgRNA_id`** — a unique name for each sgRNA (no spaces; underscores are fine)
- **`sgRNA_seq`** — the 20-nucleotide sgRNA target sequence, **without** the PAM and **without** any cloning overhangs

Additional columns (e.g. target gene, genomic coordinates, on-target score) are allowed and will be ignored by the app.

See [`example_inputs/sgRNA_seqs/`](example_inputs/sgRNA_seqs/) for ready-to-use example files in all three formats.

> **Note:** The default oligo design targets ***Sp*Cas9** sgRNAs using the ***Bsm*BI/Golden Gate** cloning system (e.g. ipUSEPR, lentiCRISPR v2, LentiCRISPRv2GFP, etc). For other cloning systems with different overhangs, change the oligo prefix/suffix fields under **Advanced options** as described in the [**Using with other vectors**](#using-with-other-vectors-and-cloning-systems) section below.

> **Tip:** The app automatically removes leading, trailing, and internal whitespace from sgRNA sequences and converts them to uppercase. Sequences copied from publications (which may include spaces every 10 bases) or pasted from PDFs are handled correctly.

---

## Step 3 — Run the app

1. Upload your sgRNA table (Step 1 in the sidebar).
2. Optional, for plasmid maps: upload your empty vector (`.gb`) and give the stuffer sequence, either typed in or as a FASTA file.
3. Click **Run**, check the preview, then download the IDT order form and the combined GenBank file.

Ready-to-use inputs for three vectors are in [`example_inputs/`](example_inputs/):

| Vector | sgRNA table | Vector file | Stuffer |
|---|---|---|---|
| ipUSEPR | `example_sgRNA_seqs.xlsx` | `vector_seqs/ipUSEPR.gb` | `GAGACGTGACCGTCTCT` (paste directly) |
| lentiCRISPR v2 ([Addgene #52961](https://www.addgene.org/52961/)) | `example_sgRNA_seqs.csv` | `vector_seqs/lentiCRISPRv2.gb` | `stuffer_seqs/lentiCRISPRv2_stuffer.fa` (1,880 bp) |
| LentiCRISPRv2GFP ([Addgene #82416](https://www.addgene.org/82416/)) | `example_sgRNA_seqs.tsv` | `vector_seqs/LentiCRISPRv2GFP.gb` | `stuffer_seqs/LentiCRISPRv2GFP_stuffer.fa` |

The **ipUSEPR** (**i**mproved lentiviral **p**lasmid with **U**6-**s**gRNA-**E**F1a-**P**uroR-Tag**R**FP) vector was developed in our [lab](https://www.mayo.edu/research/core-resources/genetic-screening-engineering-core/overview) and is our most battle-tested vector
(PMIDs: [34210975](https://pubmed.ncbi.nlm.nih.gov/34210975/),
[34373586](https://pubmed.ncbi.nlm.nih.gov/34373586/),
[36563143](https://pubmed.ncbi.nlm.nih.gov/36563143/),
[37075745](https://pubmed.ncbi.nlm.nih.gov/37075745/),
[38394203](https://pubmed.ncbi.nlm.nih.gov/38394203/),
[38316881](https://pubmed.ncbi.nlm.nih.gov/38316881/),
[39716856](https://pubmed.ncbi.nlm.nih.gov/39716856/))
for individual gene perturbations and CRISPR library screens (CRISPRko/i/a/tiling). It features a modified sgRNA backbone for enhanced sgRNA stability, a human U6 promoter driving sgRNA expression, and an EF1a promoter co-driving a puromycin resistance gene and TagRFP — enabling convenient selection of transduced cells by puromycin and rapid determination of transduction efficiency by flow cytometry. It also carries an ampicillin resistance gene for plasmid propagation in *Escherichia coli*. It is used with cells stably expressing *Sp*Cas9 or its variants, such as lentiCas9-Blast ([Addgene #52962](https://www.addgene.org/52962/)), Lenti-dCas9-KRAB-blast ([Addgene #89567](https://www.addgene.org/89567/)), and lenti dCAS-VP64_Blast ([Addgene #61425](https://www.addgene.org/61425/)).

---

## Outputs explained

The app gives you two downloads:

```
example_sgRNA_seqs_IDT_orderForm.xlsx   ← Upload this to IDT to order oligos
combined_plasmids.gb                    ← Open in SnapGene or Benchling
```

The combined `.gb` holds one record per sgRNA (e.g. `ipUSEPR_sg_Luc_474`,
`ipUSEPR_sg_Ren_55`, ...). Examples of both outputs are in
[`example_outputs/`](example_outputs/).

**IDT order form** — open in Excel to review, then upload directly to
[IDT → Oligo Entry → "BULK INPUT"](https://www.idtdna.com/site/order/oligoentry).
Each sgRNA produces two oligos: a sense (`_s`) and antisense (`_a`) strand.

**GenBank plasmid files** — open in SnapGene or Benchling to verify that the
sgRNA sequence was inserted correctly before submitting for Sanger or
whole-plasmid sequencing. Each file contains all original vector annotations
plus a new `misc_feature` labelled with the sgRNA ID at the insertion site.
The stuffer placeholder annotation is automatically removed from the output.
The construct name in each LOCUS line (e.g. `ipUSEPR_sg_Luc_474`) is
`<vector name>_<sgRNA ID>`, so SnapGene and Benchling display the correct
construct name when the file is imported.

---

## A note on the extra G in the oligo sequences

You may notice that the sense oligo begins with `CACCG` rather than just `CACC`.
The final `G` is not part of your 20-mer sgRNA sequence — it is added
deliberately for the following reason:

The **human U6 promoter**, which drives sgRNA expression in most mammalian CRISPR
vectors, is an RNA Polymerase III promoter that requires a **guanosine (G) at the
+1 transcription start site** for efficient RNA synthesis
([Ran et al., *Nat Protoc*, 2013](https://pubmed.ncbi.nlm.nih.gov/24157548/)).
When the sgRNA target sequence does not naturally begin with G, an extra G is prepended
to the sense oligo to ensure robust transcription.

Importantly, **this extra G does not significantly reduce *Sp*Cas9 on-target
cleavage activity**, as mismatches at the 5' PAM-distal end of the sgRNA are
generally better tolerated than those at the 3' PAM-proximal seed region
([Hsu et al., *Nat Biotechnol*, 2013](https://pubmed.ncbi.nlm.nih.gov/23873081/)),
though tolerance is sequence-dependent.

If your sgRNA sequence already starts with G, the extra G is still added by default,
resulting in a 21-nt sgRNA with a 5'-GG- start. This is the standard convention for
*Bsm*BI-based lentiviral CRISPR vectors, as explicitly specified in the Zhang lab's
official cloning protocol for lentiCRISPR v2 and lentiGuide-Puro
([Addgene #52963](https://www.addgene.org/52963/)), where the sense oligo is always
designed as `5'-CACCG-[20 bp target]-3'` regardless of the starting nucleotide, and
confirmed by the [Broad Institute GPP cloning protocol](https://portals.broadinstitute.org/gpp/public/dir/download?dirpath=protocols/production&filename=cloning_of_oligos_for_sgRNA_shRNA_nov2019.pdf).
This is also the convention we follow in our lab.

---

## Using with other vectors and cloning systems

The default oligo design (`CACCG` / `AAAC` / `C`) is specific to *Bsm*BI/Golden
Gate-based lentiviral vectors. The app accommodates **any** vector and cloning
system through three fields under **⚙ Advanced options**:

| Field | Default | Description |
|---|---|---|
| Sense oligo prefix | `CACCG` | 5' prefix added to the sense oligo |
| Antisense oligo prefix | `AAAC` | 5' prefix added to the antisense oligo |
| Antisense oligo suffix | `C` | 3' suffix added to the antisense oligo |

These overhangs are determined by the Type IIS restriction enzyme used to
linearize your vector. To find the correct values, check the cloning protocol
on the vector's Addgene page or its accompanying publication.

**A note on other popular vectors:** the Broad Institute **pXPR** /
**lentiGuide-Puro** series (*Bsm*BI-based) all coincidentally use the same
`CACCG` / `AAAC` / `C` overhangs as our defaults — so no customization is
needed for those vectors.

For a vector that uses different overhangs, simply change the fields.
For example (hypothetical), if your vector's cloning site generates `ACCG` / `AAAC`
overhangs with no 3' suffix (the G for U6 transcription is already embedded in the
`ACCG` overhang, so no complementary `C` is needed on the antisense oligo), set the
sense prefix to `ACCG`, the antisense prefix to `AAAC`, and leave the suffix empty.

The same **Advanced options** panel lets you set the column names for the sgRNA
sequences and IDs if your input file does not use `sgRNA_seq` / `sgRNA_id`.

---

## All inputs

| Sidebar field | Purpose |
|---|---|
| Upload sgRNA table | `.xlsx`, `.csv` or `.tsv` with the sgRNA IDs and sequences |
| Upload empty vector | Vector GenBank file (`.gb` / `.gbk`); needed only for plasmid maps |
| Stuffer sequence / FASTA | Placeholder in the vector replaced by the sgRNA. A raw sequence or a FASTA upload (the FASTA wins if both are given) |
| Vector name | Prepended to each plasmid name and LOCUS name. Default: the vector filename stem (e.g. `ipUSEPR`) |
| Combined .gb filename stem | Name of the downloaded GenBank file. Default: `combined_plasmids` |
| Advanced options | Column names and the oligo prefix/suffix described above |

---

## Supported vectors

The app works with any sgRNA cloning vector as long as you provide its
GenBank file, the stuffer sequence, and the correct cloning overhangs. The
following vectors have been tested:

| Vector | Addgene / Reference | Stuffer input | Custom overhangs needed? |
|---|---|---|---|
| ipUSEPR | PMIDs: [34210975](https://pubmed.ncbi.nlm.nih.gov/34210975/), [34373586](https://pubmed.ncbi.nlm.nih.gov/34373586/), [36563143](https://pubmed.ncbi.nlm.nih.gov/36563143/), [37075745](https://pubmed.ncbi.nlm.nih.gov/37075745/), [38394203](https://pubmed.ncbi.nlm.nih.gov/38394203/), [38316881](https://pubmed.ncbi.nlm.nih.gov/38316881/), [39716856](https://pubmed.ncbi.nlm.nih.gov/39716856/) | Short string | No (defaults) |
| lentiCRISPR v2 | [Addgene #52961](https://www.addgene.org/52961/) | FASTA file | No (defaults) |
| LentiCRISPRv2GFP | [Addgene #82416](https://www.addgene.org/82416/) | FASTA file | No (defaults) |

---

## Troubleshooting

**"Stuffer not found on either strand"**
The stuffer sequence does not match anything in your vector file. Double-check
that the sequence is correct. If your stuffer is provided as a FASTA file,
verify the uploaded file. Note: stuffers that span the circular origin of a plasmid
are not currently supported.

**"Stuffer found N times"**
The stuffer sequence is not unique in your vector. Provide a longer or more
specific stuffer sequence.

**"Stuffer sequence is palindromic"**
The stuffer's reverse complement is identical to itself, so the app cannot
determine which strand it belongs to. Use a longer, asymmetric stuffer sequence.

**"Duplicate sgRNA IDs found"**
Two or more rows in your input file share the same `sgRNA_id`. Every sgRNA must
have a unique identifier.

**"Column 'sgRNA_seq' not found"**
Your input file uses different column names. Open **⚙ Advanced options** and
enter your column names for the sgRNA sequences and IDs.

**"sgRNA contains non-ACGT characters"**
The sgRNA sequence column contains characters other than A, C, G, or T. Common
causes include ambiguity codes (N, R, Y), RNA bases (U instead of T), or
accidental inclusion of spaces or punctuation. Correct the sequence in your
input file before re-running.

**"Warning: sgRNA is N nt (expected 20 nt)"**
The app detected a sequence that is not 20 nucleotides long. *Sp*Cas9
requires a 20-nt guide for optimal activity. If the sequence was exported
from a design tool, verify that it does not include the PAM or any cloning
overhangs. The app proceeds anyway — this is a warning, not an error.

**"Warning: ends in NGG — PAM may be included"**
A sequence longer than 20 nt ends in NGG, which strongly suggests the PAM
sequence was accidentally included. Trim the sequence to 20 nt by removing
the trailing NGG before re-running.

---

## Development

The app is a Shiny for Python application in [`PyShinyApp/`](PyShinyApp/), with the core
logic in `PyShinyApp/sgrna_oligo_designer.py`. Dependencies are managed with
[uv](https://docs.astral.sh/uv/) (Python 3.12).

```bash
git clone https://github.com/anthonyknc/sgrna-oligo-designer.git
cd sgrna-oligo-designer
uv sync                                   # create .venv from uv.lock
uv run shiny run PyShinyApp/app.py --reload   # run locally
uv run pytest                             # run the tests
```

To deploy to shinyapps.io, regenerate the pinned requirements and deploy
(requires `rsconnect add` credentials):

```bash
uv export --no-dev --no-hashes --no-emit-project -o PyShinyApp/requirements.txt
uv run rsconnect deploy shiny PyShinyApp --name <account> --title sgrna-oligo-designer
```

---

## Citations

If you use this tool in your research, please cite it as:

> Chan, A.K.N. (2026). sgRNA Oligo Designer. GitHub: https://github.com/anthonyknc/sgrna-oligo-designer.

If you use our in-house vector **ipUSEPR** for sgRNA cloning, please cite:
> Yang, L.<sup>#</sup>, Chan, A.K.N.<sup>#</sup>, Miyashita, K.<sup>#</sup>, Delaney, C.D., Wang, X., Li, H., Pokharel, S.P., Li, S., Li, M., Xu, X., et al. (2021). High-resolution characterization of gene function using single-cell CRISPR tiling screen. *Nat Commun* **12**, 4063. https://doi.org/10.1038/s41467-021-24324-0. (<sup>#</sup>Co-first authors)

If you use ipUSEPR for cloning dual sgRNAs using [STAgR](https://pubmed.ncbi.nlm.nih.gov/29702666/) ([Addgene #102844](https://www.addgene.org/102844/)), please also cite:
> Chan, A.K.N., Han, L., Delaney, C.D., Wang, X., Mukhaleva, E., Li, M., Yang, L., Pokharel, S.P., Mattson, N., Garcia, M., et al. (2024). Therapeutic targeting of Tudor domains in leukemia via CRISPR-Scan Assisted Drug Discovery. *Sci Adv* **10**, eadk3127. https://doi.org/10.1126/sciadv.adk3127.

---

## License

MIT License. See [LICENSE](LICENSE) for details.

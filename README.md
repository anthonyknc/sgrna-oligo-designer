# sgRNA Oligo Designer

A command-line tool for designing sgRNA cloning oligos for CRISPR experiments.
Given a list of sgRNA sequences, it generates:

1. **An IDT bulk order form** (`.xlsx`) — ready to upload directly to [IDT](https://www.idtdna.com) for oligo synthesis
2. **Annotated plasmid maps** (`.gb`) — one GenBank file per sgRNA construct, ready to open in [SnapGene](https://www.snapgene.com) or [Benchling](https://www.benchling.com) for sequencing verification

---

## What you need before starting

- **Python 3.12 or later** — [download here](https://www.python.org/downloads/)
- **pixi** (recommended) — a fast, modern, and reproducible package manager. [Install here](https://pixi.sh/latest/)
- Or, if you prefer pip: `pip3 install openpyxl biopython`

---

## Installation

```bash
# 1. Download or clone this repository
git clone https://github.com/anthonyknc/sgrna-oligo-designer.git
cd sgrna-oligo-designer

# 2. Set up the environment (one time only)
pixi install

# 3. Activate the environment
pixi shell
```

You are now ready to run the script.

---

## Step 1 — Design your sgRNA sequences with CRISPRpick

Before using this tool, you need a list of 20-nucleotide sgRNA sequences targeting your gene(s) of interest. We recommend **CRISPRpick** from the Broad Institute:

> 🔗 https://portals.broadinstitute.org/gppx/crispick/public

CRISPRpick is free to use and requires no account. Simply enter your gene name or transcript ID, select your genome, Cas9 system (e.g. *Sp*Cas9), and the tool will return ranked sgRNA sequences with on-target efficiency scores. Select your top candidates (typically 3–6 per gene) and export them to a spreadsheet. Subsequently, the spreadsheet is tidied into the formats similar to [example_inputs/sgRNA_seqs](example_inputs/sgRNA_seqs).

---

## Step 2 — Prepare your input file

Prepare a spreadsheet (Excel `.xlsx`, `.csv`, or `.tsv`) with at least these two columns:

| sgRNA_id | sgRNA_seq |
|---|---|
| sg_Luc_474 | GATTCTAAAACGGATTACCA |
| sg_tRFP657_294 | GTCACCACATACGAAGACGG |

- **`sgRNA_id`** — a unique name for each sgRNA (no spaces; underscores are fine)
- **`sgRNA_seq`** — the 20-nucleotide sgRNA sequence, **without** the PAM and **without** any prefixes

Additional columns (e.g. target gene, genomic coordinates, on-target score) are allowed and will be ignored by the script.

See [`example_inputs/sgRNA_seqs/`](example_inputs/sgRNA_seqs/) for ready-to-use example files in all three formats.

> **Note:** This script is designed for ***Sp*Cas9** sgRNAs using the ***Bsm*BI/Golden Gate** cloning system (e.g. ipUSEPR, lentiCRISPRv2). For other nucleases or cloning systems, adjust the `--top-prefix`, `--bot-prefix`, and `--bot-suffix` options accordingly.

---

## Step 3 — Run the script

### Example 1 — ipUSEPR vector (Excel input)

The **ipUSEPR** (**i**mproved lentiviral **p**lasmid with **U**6-**s**gRNA-**E**F1a-**P**uroR-Tag**R**FP) vector was developed in our [lab](https://www.mayo.edu/research/core-resources/genetic-screening-engineering-core/overview) and is our most battle-tested vector for individual gene perturbations and CRISPR library screens (CRISPRko/i/a/tiling). It features a modified sgRNA backbone for enhanced sgRNA stability, a U6 promoter driving sgRNA expression, and an EF1a promoter driving co-expression of a puromycin resistance gene and TagRFP — enabling dual selection of transduced cells by antibiotic resistance and red fluorescence. It also carries an ampicillin resistance gene for plasmid amplification in *Escherichia coli*. It is used with cells stably expressing *Sp*Cas9 or its variants, such as lentiCas9-Blast ([Addgene #52962](https://www.addgene.org/52962/)), Lenti-dCas9-KRAB-blast ([Addgene #89567](https://www.addgene.org/89567/)), and lenti-dCAS-VP64-Blast ([Addgene #61425](https://www.addgene.org/61425/))
(PMIDs: [34210975](https://pubmed.ncbi.nlm.nih.gov/34210975/),
[34373586](https://pubmed.ncbi.nlm.nih.gov/34373586/),
[36563143](https://pubmed.ncbi.nlm.nih.gov/36563143/),
[37075745](https://pubmed.ncbi.nlm.nih.gov/37075745/),
[38394203](https://pubmed.ncbi.nlm.nih.gov/38394203/),
[38316881](https://pubmed.ncbi.nlm.nih.gov/38316881/),
[39716856](https://pubmed.ncbi.nlm.nih.gov/39716856/)).
The stuffer sequence is short enough to paste directly as a string.

```bash
python3 sgrna_oligo_designer.py example_inputs/sgRNA_seqs/example_sgRNA_seqs.xlsx \
    --vector example_inputs/vector_seqs/ipUSEPR.gb \
    --stuffer GAGACGTGACCGTCTCT \
    --output example_outputs/ \
    --gb-output example_outputs/plasmids/
```

### Example 2 — lentiCRISPRv2 vector (CSV input)

**lentiCRISPRv2** ([Addgene #52961](https://www.addgene.org/52961/)) is one of the most widely used all-in-one lentiviral CRISPR vectors. Its stuffer sequence is long (~1,700 bp), so it is provided as a FASTA file — the script detects this automatically.

```bash
python3 sgrna_oligo_designer.py example_inputs/sgRNA_seqs/example_sgRNA_seqs.csv \
    --vector example_inputs/vector_seqs/lentiCRISPRv2.gb \
    --stuffer example_inputs/stuffer_seqs/lentiCRISPRv2_stuffer.fa \
    --output example_outputs/ \
    --gb-output example_outputs/plasmids/
```

### Example 3 — LentiCRISPRv2GFP vector (TSV input)

**LentiCRISPRv2GFP** ([Addgene #82416](https://www.addgene.org/82416/)) is a lentiCRISPRv2 variant with a GFP reporter, useful for tracking transduced cells by flow cytometry and growth competition assays. It also uses a FASTA stuffer file.

```bash
python3 sgrna_oligo_designer.py example_inputs/sgRNA_seqs/example_sgRNA_seqs.tsv \
    --vector example_inputs/vector_seqs/LentiCRISPRv2GFP.gb \
    --stuffer example_inputs/stuffer_seqs/LentiCRISPRv2GFP_stuffer.fa \
    --output example_outputs/ \
    --gb-output example_outputs/plasmids/
```

---

## Outputs explained

After running the script, you will find:

```
example_outputs/
├── example_sgRNA_seqs_IDT_orderForm.xlsx   ← Upload this to IDT to order oligos
└── plasmids/
    ├── ipUSEPR_sg_Luc_474.gb               ← Open in SnapGene or Benchling
    ├── ipUSEPR_sg_Ren_55.gb
    └── ...
```

**IDT order form** — open in Excel to review, then upload directly to
[IDT → "Order" → "Bulk Input"](https://www.idtdna.com/site/order/oligoentry).
Each sgRNA produces two oligos: a sense (`_s`) and antisense (`_a`) strand.

**GenBank plasmid files** — open in SnapGene or Benchling to verify that the
sgRNA sequence was inserted correctly before submitting for Sanger or
whole-plasmid sequencing. Each file contains all original vector annotations
plus a new `misc_feature` labelled with the sgRNA ID at the insertion site.

---

## A note on the extra G in the oligo sequences

You may notice that the sense oligo begins with `CACCG` rather than just `CACC`.
The final `G` is not part of your 20-mer sgRNA sequence — it is added
deliberately for the following reason:

The **human U6 promoter**, which drives sgRNA expression in most mammalian CRISPR
vectors, is an RNA Polymerase III promoter that requires a **guanosine (G) at the
+1 transcription start site** for efficient RNA synthesis
([Ran et al., *Nat Protoc*, 2013](https://pubmed.ncbi.nlm.nih.gov/24157548/)).
When the sgRNA sequence does not naturally begin with G, an extra G is prepended
to ensure robust transcription.

Importantly, **this extra G does not significantly affect *Sp*Cas9 on-target
cleavage activity** for wild-type *Sp*Cas9, as position 1 mismatches at the 5'
end of the sgRNA are well tolerated
([Hsu et al., *Nat Biotechnol*, 2013](https://pubmed.ncbi.nlm.nih.gov/23873081/);
[Ran et al., *Nat Protoc*, 2013](https://pubmed.ncbi.nlm.nih.gov/24157548/)).
If your sgRNA sequence already starts with G, the extra G is still added by
default (as the Zhang lab recommends), resulting in a 5'-GG- start — this has
been shown to have minimal impact on activity and may marginally improve
specificity ([Doench et al., *Nat Biotechnol*, 2016](https://pubmed.ncbi.nlm.nih.gov/26780180/)).

> If you are using a **high-fidelity Cas9 variant** (e.g. eSpCas9, SpCas9-HF1),
> note that these variants are more sensitive to 5' mismatches. In that case,
> preferentially select sgRNA sequences that naturally begin with G to avoid
> the mismatch altogether.

---

## All options

```
usage: sgrna_oligo_designer [-h] [--output OUTPUT]
                            [--vector VECTOR] [--stuffer STUFFER]
                            [--multi-gb] [--gb-output GB_OUTPUT]
                            [--vector-name VECTOR_NAME]
                            [--seq-col SEQ_COL] [--id-col ID_COL]
                            [--sense-prefix SENSE_PREFIX]
                            [--antisense-prefix ANTISENSE_PREFIX]
                            [--antisense-suffix ANTISENSE_SUFFIX]
                            input

positional arguments:
  input                     Input file (.xlsx, .csv, or .tsv)

optional arguments:
  --output, -o              Output directory or path for the IDT order form
                            Default: same folder as the input file
  --vector, -v              Empty vector GenBank file (.gb). Required for plasmid output
  --stuffer                 Stuffer sequence to replace with the sgRNA 20-mer.
                            Accepts a raw sequence string or a FASTA file path
                            (auto-detected). Default: 20 N's (NNNNNNNNNNNNNNNNNNNN)
  --multi-gb                Combine all plasmid records into a single .gb file
                            instead of one file per sgRNA
  --gb-output               Output directory for .gb plasmid files
                            Default: same folder as the input file
  --vector-name             Custom name prepended to each plasmid filename.
                            Default: stem of the --vector filename (e.g. "ipUSEPR")
  --seq-col                 Column name for sgRNA sequences (default: sgRNA_seq)
  --id-col                  Column name for sgRNA identifiers (default: sgRNA_id)
  --sense-prefix            5' prefix for the sense oligo (default: CACCG)
  --antisense-prefix        5' prefix for the antisense oligo (default: AAAC)
  --antisense-suffix        3' suffix for the antisense oligo (default: C)
```

---

## Supported vectors

The script works with **any** *Bsm*BI/Golden Gate-compatible vector as long as you
provide its GenBank file and the stuffer sequence to replace. The following
vectors have been tested:

| Vector | Addgene / Reference | Stuffer input |
|---|---|---|
| ipUSEPR | PMIDs: 34210975, 34373586, 36563143, 37075745, 38394203, 38316881, 39716856 | Short string |
| lentiCRISPRv2 | [Addgene #52961](https://www.addgene.org/52961/) | FASTA file |
| LentiCRISPRv2GFP | [Addgene #82416](https://www.addgene.org/82416/) | FASTA file |

---

## Troubleshooting

**"Stuffer not found on either strand"**
The stuffer sequence does not match anything in your vector file. Double-check
that the sequence is correct. If your stuffer is provided as a FASTA file,
verify the file path. Note: stuffers that span the circular origin of a plasmid
are not currently supported.

**"Stuffer found N times"**
The stuffer sequence is not unique in your vector. Provide a longer or more
specific stuffer sequence.

**"Stuffer sequence is palindromic"**
The stuffer's reverse complement is identical to itself, so the script cannot
determine which strand it belongs to. Use a longer, asymmetric stuffer sequence.

**"Duplicate sgRNA IDs found"**
Two or more rows in your input file share the same `sgRNA_id`. Every sgRNA must
have a unique identifier.

**"Column 'sgRNA_seq' not found"**
Your input file uses different column names. Use `--seq-col` and `--id-col` to
specify the correct column names, e.g.:
```bash
python3 sgrna_oligo_designer.py input.xlsx --seq-col guide_sequence --id-col guide_id
```

---

## Citations

If you use this tool in your research, please cite it as:

> Chan, A.K.N. (2026). sgRNA Oligo Designer. GitHub: https://github.com/anthonyknc/sgrna-oligo-designer.

If you use our in-house plasmid, ipUSEPR for sgRNA cloning, please cite it as:
> Yang, L.<sup>#</sup>, Chan, A.K.N.<sup>#</sup>, Miyashita, K.<sup>#</sup>, Delaney, C.D., Wang, X., Li, H., Pokharel, S.P., Li, S., Li, M., Xu, X., et al. (2021). High-resolution characterization of gene function using single-cell CRISPR tiling screen. Nat Commun *12*, 4063. 10.1038/s41467-021-24324-0. (<sup>#</sup>Co-first authors)

If you use ipUSEPR for duo sgRNAs cloning using [STAgR](https://pubmed.ncbi.nlm.nih.gov/29702666/) ([Addgene #102844](https://www.addgene.org/102844/)), please cite it as:
> Chan, A.K.N., Han, L., Delaney, C.D., Wang, X., Mukhaleva, E., Li, M., Yang, L., Pokharel, S.P., Mattson, N., Garcia, M., et al. (2024). Therapeutic targeting Tudor domains in leukemia via CRISPR-Scan Assisted Drug Discovery. Sci Adv *10*, eadk3127. 10.1126/sciadv.adk3127.

---

## License

MIT License. See [LICENSE](LICENSE) for details.

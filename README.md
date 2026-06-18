# sgRNA Oligo Designer

A command-line tool for designing sgRNA cloning oligos for CRISPR experiments.
Given a list of sgRNA sequences, it generates:

1. **An IDT bulk order form** (`.xlsx`) — ready to upload directly to [IDT](https://www.idtdna.com) for oligo synthesis
2. **Annotated plasmid maps** (`.gb`) — one GenBank file per sgRNA construct, ready to open in [SnapGene](https://www.snapgene.com) or [Benchling](https://www.benchling.com) for sequencing verification

---

## What you need before starting

- **Python 3.12 or later** — [download here](https://www.python.org/downloads/)
- **pixi** (recommended) — a fast, easy package manager for scientific Python. [Install instructions here](https://pixi.sh/latest/)
- Or, if you prefer: `pip install openpyxl biopython`

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

## Input file format

Prepare a spreadsheet (Excel `.xlsx`, `.csv`, or `.tsv`) with at least these two columns:

| sgRNA_id | sgRNA_seq |
|---|---|
| sg_Luc_474 | GATTCTAAAACGGATTACCA |
| sg_tRFP657_294 | GTCACCACATACGAAGACGG |

- **`sgRNA_id`** — a unique name for each sgRNA (no spaces; underscores are fine)
- **`sgRNA_seq`** — the 20-nucleotide sgRNA sequence, **without** the PAM

Additional columns (e.g. target gene, genomic coordinates, on-target score) are allowed and will be ignored by the script.

See [`example_inputs/sgRNA_seqs/`](example_inputs/sgRNA_seqs/) for ready-to-use example files in all three formats.

---

## Quick start examples

The following examples use the provided example input files and reproduce the
outputs in [`example_outputs/`](example_outputs/).

### Example 1 — ipUSEPR vector (Excel input)

The **ipUSEPR** vector is most battle-testd sgRNA expression system in our [lab](https://www.mayo.edu/research/core-resources/genetic-screening-engineering-core/overview) for individual gene pertubations and CRISPR library screens using CRISPRko/i/a/tiling. It is used with cells engineered with SpCas9 and its variants, e.g. lentiCas9-Blast [Addgene #52962](https://www.addgene.org/52962/), Lenti-dCas9-KRAB-blast [Addgene #89567](https://www.addgene.org/89567/), and lenti dCAS-VP64_Blast [Addgene #61425](https://www.addgene.org/61425/), etc.
(PMIDs: [34210975](https://pubmed.ncbi.nlm.nih.gov/34210975/),
[34373586](https://pubmed.ncbi.nlm.nih.gov/34373586/),
[36563143](https://pubmed.ncbi.nlm.nih.gov/36563143/),
[37075745](https://pubmed.ncbi.nlm.nih.gov/37075745/),
[38394203](https://pubmed.ncbi.nlm.nih.gov/38394203/),
[38316881](https://pubmed.ncbi.nlm.nih.gov/38316881/),
[39716856](https://pubmed.ncbi.nlm.nih.gov/39716856/)).
The stuffer is short enough to paste directly as a string.

```bash
python3 sgrna_oligo_designer.py example_inputs/sgRNA_seqs/example_sgRNA_seqs.xlsx \
    --vector example_inputs/vector_seqs/ipUSEPR.gb \
    --stuffer GAGACGTGACCGTCTCT \
    --output example_outputs/ \
    --gb-output example_outputs/plasmids/
```

### Example 2 — lentiCRISPRv2 vector (CSV input)

**lentiCRISPRv2** ([Addgene #52961](https://www.addgene.org/52961/)) is one of
the most widely used all-in-one lentiviral CRISPR vectors. Its stuffer sequence
is long (~1,700 bp), so it is provided as a FASTA file — the script detects
this automatically.

```bash
python3 sgrna_oligo_designer.py example_inputs/sgRNA_seqs/example_sgRNA_seqs.csv \
    --vector example_inputs/vector_seqs/lentiCRISPRv2.gb \
    --stuffer example_inputs/stuffer_seqs/lentiCRISPRv2_stuffer.fa \
    --output example_outputs/ \
    --gb-output example_outputs/plasmids/
```

### Example 3 — LentiCRISPRv2GFP vector (TSV input)

**LentiCRISPRv2GFP** ([Addgene #82416](https://www.addgene.org/82416/)) is a
lentiCRISPRv2 variant with a GFP reporter, useful for tracking transduced cells
by flow cytometry and growth competition assay. Also uses a FASTA stuffer file.

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
[IDT's order page –> "BULK INPUT"](https://www.idtdna.com/site/order/oligoentry).
Each sgRNA produces two oligos: a sense (`_s`) and antisense (`_a`) strand.

**GenBank plasmid files** — open in SnapGene or Benchling to verify that the
sgRNA was inserted correctly before sending for Sanger or whole-plasmid
sequencing. Each file contains all original vector annotations plus a new
`misc_feature` labelled with the sgRNA ID at the insertion site.

---

## All options

```
usage: sgrna_oligo_designer [-h] [--output OUTPUT]
                            [--vector VECTOR] [--stuffer STUFFER]
                            [--multi-gb] [--gb-output GB_OUTPUT]
                            [--vector-name VECTOR_NAME]
                            [--seq-col SEQ_COL] [--id-col ID_COL]
                            [--top-prefix TOP_PREFIX]
                            [--bot-prefix BOT_PREFIX]
                            [--bot-suffix BOT_SUFFIX]
                            input

positional arguments:
  input                     Input file (.xlsx, .csv, or .tsv)

optional arguments:
  --output, -o              Output directory or path for IDT order form Excel file
                            Default: same folder as input file
  --vector, -v              Empty vector GenBank file (.gb). Required for plasmid output
  --stuffer                 Stuffer sequence to replace with sgRNA. Accepts a raw
                            sequence string or a path to a FASTA file (auto-detected)
                            Default: 20 N's (NNNNNNNNNNNNNNNNNNNN)
  --multi-gb                Combine all plasmid records into a single .gb file
                            instead of one file per sgRNA
  --gb-output               Output directory for .gb plasmid files
                            Default: same folder as input file
  --vector-name             Custom name appended to each plasmid filename
                            Default: stem of the --vector filename (e.g. "ipUSEPR")
  --seq-col                 Column name for sgRNA sequences (default: sgRNA_seq)
  --id-col                  Column name for sgRNA identifiers (default: sgRNA_id)
  --top-prefix              5' prefix for sense oligo (default: CACCG)
  --bot-prefix              5' prefix for antisense oligo (default: AAAC)
  --bot-suffix              3' suffix for antisense oligo (default: C)
```

> **Default oligo design targets BsmBI-based vectors** (e.g. ipUSEPR, lentiCRISPRv2).
> For other cloning systems, adjust `--top-prefix`, `--bot-prefix`, and `--bot-suffix`
> accordingly.

---

## Supported vectors

The script works with **any** vector as long as you provide its GenBank file
and the stuffer sequence to replace. The following vectors have been tested:

| Vector | Addgene / Reference | Stuffer |
|---|---|---|
| ipUSEPR | PMIDs: 34210975, 34373586, 36563143, 37075745, 38394203, 38316881, 39716856 | Short string |
| lentiCRISPRv2 | [Addgene #52961](https://www.addgene.org/52961/) | FASTA file |
| LentiCRISPRv2GFP | [Addgene #82416](https://www.addgene.org/82416/) | FASTA file |

---

## Troubleshooting

**"Stuffer not found on either strand"**
The stuffer sequence does not match anything in your vector file. Double-check
that the sequence is correct, or try the reverse complement. If your stuffer is
in a FASTA file, make sure the file path is correct.

**"Stuffer found N times"**
The stuffer sequence is not unique in your vector. Provide a longer or more
specific stuffer sequence.

**"Duplicate sgRNA IDs found"**
Two or more rows in your input file share the same `sgRNA_id`. Each sgRNA must
have a unique name.

**"Column 'sgRNA_seq' not found"**
Your input file uses different column names. Use `--seq-col` and `--id-col` to
specify the correct column names.

---

## Citation

If you use this tool in your research, please cite it as:

> Chan, A.K.N. (2026). sgRNA Oligo Designer.
> GitHub: https://github.com/anthonyknc/sgrna-oligo-designer

---

## License

MIT License. See [LICENSE](LICENSE) for details.

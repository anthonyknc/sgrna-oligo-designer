# Getting started

## Step 1 — Design your sgRNA sequences with CRISPick

Before using this tool, you need a list of 20-nucleotide sgRNA target sequences for your gene(s) of interest. We recommend **CRISPick** from Broad Institute:

> 🔗 https://portals.broadinstitute.org/gppx/crispick/public

CRISPick is free to use and requires no account. Simply enter your gene name or transcript ID, select your genome and Cas9 system (e.g. *Sp*Cas9), and the tool will return ranked sgRNA sequences with on-target efficiency scores. Select your top candidates (typically 3–6 per gene) and export them to a spreadsheet, then reformat the exported data to match the input format shown in the examples at [`example_inputs/sgRNA_seqs/`](https://github.com/anthonyknc/sgrna-oligo-designer/tree/main/example_inputs/sgRNA_seqs).


## Step 2 — Prepare your input file

Prepare a spreadsheet (Excel `.xlsx`, `.csv`, or `.tsv`) with at least these two columns:

| sgRNA_id | sgRNA_seq |
|---|---|
| sg_Luc_474 | GATTCTAAAACGGATTACCA |
| sg_tRFP657_294 | GTCACCACATACGAAGACGG |

- **`sgRNA_id`** — a unique name for each sgRNA (no spaces; underscores are fine)
- **`sgRNA_seq`** — the 20-nucleotide sgRNA target sequence, **without** the PAM and **without** any cloning overhangs

Additional columns (e.g. target gene, genomic coordinates, on-target score) are allowed and will be ignored by the app.

See [`example_inputs/sgRNA_seqs/`](https://github.com/anthonyknc/sgrna-oligo-designer/tree/main/example_inputs/sgRNA_seqs) for ready-to-use example files in all three formats.

> **Note:** The tool designs oligos for ***Sp*Cas9** sgRNAs cloned by Type IIS/Golden Gate digestion (*Bsm*BI, Esp3I, *Bsa*I or *Bbs*I). You never type overhangs: they are read from the vector, as explained under **How it works**.

> **Tip:** The app automatically removes leading, trailing, and internal whitespace from sgRNA sequences and converts them to uppercase. Sequences copied from publications (which may include spaces every 10 bases) or pasted from PDFs are handled correctly.


## Step 3 — Run the app

1. Upload your sgRNA table (Step 1 in the sidebar).
2. Pick your **vector** (ipUSEPR, lentiCRISPR v2, LentiCRISPRv2GFP, or upload your own `.gb`) and the **restriction enzyme** used to clone into it. The app immediately shows the cloning site it found: the two overhangs and how many bases are replaced.
3. Click **Run**, check the preview and any warnings, then download the IDT order form and the combined GenBank file.

Ready-to-use sgRNA tables are in [`example_inputs/sgRNA_seqs/`](https://github.com/anthonyknc/sgrna-oligo-designer/tree/main/example_inputs/sgRNA_seqs) in three formats, and the three vectors are bundled with the app (and in [`PyShinyApp/vectors/`](https://github.com/anthonyknc/sgrna-oligo-designer/tree/main/PyShinyApp/vectors)). All three are *Bsm*BI vectors.

The **ipUSEPR** (**i**mproved lentiviral **p**lasmid with **U**6-**s**gRNA-**E**F1a-**P**uroR-Tag**R**FP) vector was developed in our [lab](https://www.mayo.edu/research/core-resources/genetic-screening-engineering-core/overview) and is our most battle-tested vector
(PMIDs: [34210975](https://pubmed.ncbi.nlm.nih.gov/34210975/),
[34373586](https://pubmed.ncbi.nlm.nih.gov/34373586/),
[36563143](https://pubmed.ncbi.nlm.nih.gov/36563143/),
[37075745](https://pubmed.ncbi.nlm.nih.gov/37075745/),
[38394203](https://pubmed.ncbi.nlm.nih.gov/38394203/),
[38316881](https://pubmed.ncbi.nlm.nih.gov/38316881/),
[39716856](https://pubmed.ncbi.nlm.nih.gov/39716856/))
for individual gene perturbations and CRISPR library screens (CRISPRko/i/a/tiling). It features a modified sgRNA backbone for enhanced sgRNA stability, a human U6 promoter driving sgRNA expression, and an EF1a promoter co-driving a puromycin resistance gene and TagRFP — enabling convenient selection of transduced cells by puromycin and rapid determination of transduction efficiency by flow cytometry. It also carries an ampicillin resistance gene for plasmid propagation in *Escherichia coli*. It is used with cells stably expressing *Sp*Cas9 or its variants, such as lentiCas9-Blast ([Addgene #52962](https://www.addgene.org/52962/)), Lenti-dCas9-KRAB-blast ([Addgene #89567](https://www.addgene.org/89567/)), and lenti dCAS-VP64_Blast ([Addgene #61425](https://www.addgene.org/61425/)).


## All inputs

| Sidebar field | Purpose |
|---|---|
| Upload sgRNA table | `.xlsx`, `.csv` or `.tsv` with the sgRNA IDs and sequences |
| Vector | A bundled vector, or upload your own circular GenBank file (`.gb` / `.gbk`) |
| Restriction enzyme | *Bsm*BI, Esp3I, *Bsa*I or *Bbs*I; sets the cuts, overhangs and stuffer |
| Vector name | Prepended to each plasmid name and LOCUS name. Default: the vector name or uploaded file name |
| Combined .gb filename stem | Name of the downloaded GenBank file. Default: `combined_plasmids` |
| Advanced options | Column names; whether to add the 5' G (always / only if the sgRNA doesn't start with G / never); sgRNA orientation |


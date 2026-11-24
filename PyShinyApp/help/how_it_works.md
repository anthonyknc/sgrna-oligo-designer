# How it works

## How the cloning is simulated

Instead of asking you for overhang sequences, the app digests your vector in silico:

1. It finds the enzyme's two recognition sites (on either strand) and works out where each one cuts. Type IIS enzymes cut outside their site: *Bsm*BI and *Bsa*I cut `N1/N5` downstream of the site, *Bbs*I cuts `N2/N6`.
2. The 4-nt overhangs left on the vector at those cuts are the overhangs your oligos need.
3. The stuffer is the shorter piece between the two cuts. It is replaced by the sgRNA, so you no longer supply a stuffer sequence.
4. Oligos are built as `sense = left overhang + [G] + sgRNA` and `antisense = reverse complement of ([G] + sgRNA + right overhang)`.

For ipUSEPR with *Bsm*BI the overhangs come out as `CACC` and `GTTT`, which reproduces the familiar `CACCG-sgRNA` / `AAAC-rc(sgRNA)-C` oligos.

Both common vector layouts work: recognition sites inside the stuffer, lost on cutting (lentiCRISPR, ipUSEPR), and sites that stay in the backbone. The sgRNA orientation (which strand of the GenBank file is transcribed) is auto-detected from the Cas9 scaffold start (`GTTT`), and you can override it under **⚙ Advanced options**.

The app checks the following and tells you if one fails:
- The vector must be circular and have **exactly two** sites for the chosen enzyme. If you pick the wrong enzyme, the message lists the enzymes that do cut your vector exactly twice.
- Each sgRNA is checked for the enzyme's recognition site, on either strand and including sites created at the junctions. The enzyme would cut such a guide during cloning.

**Supported enzymes:** *Bsm*BI, Esp3I, *Bsa*I and *Bbs*I. (*Bsa*I does not cut the bundled vectors exactly twice, so choose it only with a *Bsa*I vector of your own.)


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

If your sgRNA sequence already starts with G, the extra G is still added by default (you can change this under **Advanced options**),
resulting in a 21-nt sgRNA with a 5'-GG- start. This is the standard convention for
*Bsm*BI-based lentiviral CRISPR vectors, as explicitly specified in the Zhang lab's
official cloning protocol for lentiCRISPR v2 and lentiGuide-Puro
([Addgene #52963](https://www.addgene.org/52963/)), where the sense oligo is always
designed as `5'-CACCG-[20 bp target]-3'` regardless of the starting nucleotide, and
confirmed by the [Broad Institute GPP cloning protocol](https://portals.broadinstitute.org/gpp/public/dir/download?dirpath=protocols/production&filename=cloning_of_oligos_for_sgRNA_shRNA_nov2019.pdf).
This is also the convention we follow in our lab.


## Outputs explained

The app gives you two downloads:

```
example_sgRNA_seqs_IDT_orderForm.xlsx   ← Upload this to IDT to order oligos
combined_plasmids.gb                    ← Open in SnapGene or Benchling
```

The combined `.gb` holds one record per sgRNA (e.g. `ipUSEPR_sg_Luc_474`,
`ipUSEPR_sg_Ren_55`, ...). Examples of both outputs are in
[`example_outputs/`](https://github.com/anthonyknc/sgrna-oligo-designer/tree/main/example_outputs).

**IDT order form** — open in Excel to review, then upload directly to
[IDT → Oligo Entry → "BULK INPUT"](https://www.idtdna.com/site/order/oligoentry).
Each sgRNA produces two oligos: a sense (`_s`) and antisense (`_a`) strand.

**GenBank plasmid files** — open in SnapGene or Benchling to verify that the
sgRNA sequence was inserted correctly before submitting for Sanger or
whole-plasmid sequencing. Each file contains all original vector annotations
plus a new `misc_feature` labelled with the sgRNA ID at the insertion site.
Annotations on the replaced stuffer are removed, and annotations that overlap it are trimmed.
The construct name in each LOCUS line (e.g. `ipUSEPR_sg_Luc_474`) is
`<vector name>_<sgRNA ID>`, so SnapGene and Benchling display the correct
construct name when the file is imported.


## Supported vectors

The app works with any circular sgRNA cloning vector that has exactly two sites for one of the supported enzymes. The following are tested (the example outputs in [`example_outputs/`](https://github.com/anthonyknc/sgrna-oligo-designer/tree/main/example_outputs) are generated from them):

| Vector | Addgene / Reference | Enzyme | sgRNA orientation in the file |
|---|---|---|---|
| ipUSEPR | PMIDs: [34210975](https://pubmed.ncbi.nlm.nih.gov/34210975/), [34373586](https://pubmed.ncbi.nlm.nih.gov/34373586/), [36563143](https://pubmed.ncbi.nlm.nih.gov/36563143/), [37075745](https://pubmed.ncbi.nlm.nih.gov/37075745/), [38394203](https://pubmed.ncbi.nlm.nih.gov/38394203/), [38316881](https://pubmed.ncbi.nlm.nih.gov/38316881/), [39716856](https://pubmed.ncbi.nlm.nih.gov/39716856/) | *Bsm*BI | forward |
| lentiCRISPR v2 | [Addgene #52961](https://www.addgene.org/52961/) | *Bsm*BI | forward |
| LentiCRISPRv2GFP | [Addgene #82416](https://www.addgene.org/82416/) | *Bsm*BI | reverse |


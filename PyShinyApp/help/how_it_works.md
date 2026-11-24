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

**Supported enzymes:** *Bsm*BI/Esp3I, *Bbs*I/BpiI, *Bsa*I and SapI are listed first (BsmBI and Esp3I, and BbsI and BpiI, are the same enzyme under different names). Choose **Other…** to search every commercially available Type IIS enzyme in Biopython's restriction enzyme database that recognises a plain A/C/G/T site, cuts outside it, and leaves a 3 to 5 nt 5' overhang (for example BspQI, LguI, BtgZI, AarI or FokI). Enzymes with 1 to 2 nt overhangs, blunt or two-sided cutters, and enzymes with degenerate sites are not offered in this mode.

For an enzyme from the **Other…** list, the app only checks the cut geometry. Check yourself that the enzyme suits your vector (for example methylation sensitivity) and follow the supplier's protocol. For overhangs shorter than 3 nt the sgRNA orientation can't be auto-detected, and with 3-nt overhangs the scaffold check uses only the first three bases (`GTT`).

(*Bsa*I does not cut the bundled vectors exactly twice, so choose it only with a *Bsa*I vector of your own.)


## A note on the extra G in the oligo sequences

You may notice that the sense oligo begins with `CACCG` rather than just `CACC`.
The final `G` is not part of your 20-mer sgRNA sequence. It is added on purpose.

The **human U6 promoter**, which drives sgRNA expression in most mammalian CRISPR
vectors, is an RNA polymerase III promoter that **strongly prefers a guanosine (G) as the
first transcribed base** (+1)
([Ran et al., *Nat Protoc*, 2013](https://pubmed.ncbi.nlm.nih.gov/24157548/)).
A 5' G therefore helps ensure efficient transcription.

**By default the app adds this G to every sgRNA, whether or not the sequence already
starts with G.** A sgRNA that starts with G becomes a 21-nt guide with a 5'-GG- start.
This is the standard convention for *Bsm*BI-based lentiviral CRISPR vectors, as specified in
the Zhang lab's cloning protocol for lentiCRISPR v2 and lentiGuide-Puro
([Addgene #52963](https://www.addgene.org/52963/)), where the sense oligo is always
designed as `5'-CACCG-[20 bp target]-3'`, and in the
[Broad Institute GPP cloning protocol](https://portals.broadinstitute.org/gpp/public/dir/download?dirpath=protocols/production&filename=cloning_of_oligos_for_sgRNA_shRNA_nov2019.pdf).
This is also the convention we follow in our lab.

Under **Advanced options** you can instead add the G only when the sgRNA does not already
start with one, or never add it, for example if your sgRNA is driven by a promoter other
than U6. Other RNA polymerase III promoters (such as H1 or 7SK) may prefer a different
first base, so check the requirements of yours.

The extra 5' G is widely used and is generally well tolerated in practice. This fits with
the observation that SpCas9 tolerates mismatches at the PAM-distal 5' end of the guide
better than at the PAM-proximal seed region
([Hsu et al., *Nat Biotechnol*, 2013](https://pubmed.ncbi.nlm.nih.gov/23873081/)).
That study looked at mismatches rather than an added base, so treat it as supporting
evidence rather than proof, and note that tolerance depends on the sequence.

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


## Vectors cut by two different enzymes

Some sgRNA vectors are cut by two different enzymes rather than one enzyme at two sites. A common example is BstXI + BlpI (the pMCB306 / pMCB320 / CROP-seq family); another is AgeI + EcoRI (the pLKO.1 shRNA family). Choose **Two different enzymes** under Cloning enzymes and pick both enzymes. The lists contain only the enzymes that cut your vector **exactly once**, and the order does not matter.

Any commercially available enzyme can be used here, including palindromic and degenerate sites and enzymes that leave a 3' overhang (BstXI leaves a 4-nt 3' overhang, BlpI a 3-nt 5' overhang). The app records, for each cut, which strand is cut where, and builds each oligo as the part of the insert that the vector ends do not already supply. Cloning the vector's own original guide back in therefore restores the original plasmid.

**Constant bases around the spacer.** With two enzymes the cuts do not say where the spacer goes, because the vector keeps some constant bases inside the removed piece (for pMCB306: `TTG` before the spacer and `GTTTAAGAGC` after it, the start of the scaffold). The app works these out by finding the Cas9 scaffold start in the removed piece and assuming a 20-nt spacer slot, and shows them under **Cloning site**. Under **Advanced options** you can type your own values (copy them from your vector's oligo design instructions); a blank box means "use the automatic ones". If the scaffold is not found, enter them yourself.

Two-enzyme mode has only been tested on one real vector (pMCB306) and on synthetic vectors, so compare the designed oligos with your vector's published oligo design before ordering. It currently designs sgRNAs only; shRNA hairpin vectors such as pLKO.1 need a different insert and are not supported yet.

## Do the oligos need to be phosphorylated?

No. The oligos are ordered without a 5' phosphate, and that is fine. The digested vector
ends carry 5' phosphates, so the ligase joins each vector end to the 3' end of the annealed
oligos. The two remaining nicks (at the oligos' 5' ends) are repaired in the bacteria after
transformation. Some published protocols add a T4 PNK step to phosphorylate the oligos
first; that is optional here.

## Supported vectors

The app works with any circular sgRNA cloning vector that has exactly two sites for one of the supported enzymes. The following are tested (the example outputs in [`example_outputs/`](https://github.com/anthonyknc/sgrna-oligo-designer/tree/main/example_outputs) are generated from them):

| Vector | Addgene / Reference | Enzyme | sgRNA orientation in the file |
|---|---|---|---|
| ipUSEPR | PMIDs: [34210975](https://pubmed.ncbi.nlm.nih.gov/34210975/), [34373586](https://pubmed.ncbi.nlm.nih.gov/34373586/), [36563143](https://pubmed.ncbi.nlm.nih.gov/36563143/), [37075745](https://pubmed.ncbi.nlm.nih.gov/37075745/), [38394203](https://pubmed.ncbi.nlm.nih.gov/38394203/), [38316881](https://pubmed.ncbi.nlm.nih.gov/38316881/), [39716856](https://pubmed.ncbi.nlm.nih.gov/39716856/) | *Bsm*BI | forward |
| lentiCRISPR v2 | [Addgene #52961](https://www.addgene.org/52961/) | *Bsm*BI | forward |
| LentiCRISPRv2GFP | [Addgene #82416](https://www.addgene.org/82416/) | *Bsm*BI | reverse |


# Bundled vectors

These GenBank files are the public sequence records of the plasmids below, offered as presets in the
app and used as test data. The plasmids themselves are distributed by Addgene (or the depositor) under
the depositor's MTA; only the sequence records are kept here. Cite the paper that describes the plasmid
when you use it.

| File | Plasmid | Reference | Cloning enzyme(s) |
|---|---|---|---|
| `ipUSEPR.gb` | ipUSEPR | see the Citing page of the app | BsmBI |
| `lentiCRISPRv2.gb` | lentiCRISPR v2, Addgene 52961 | Sanjana et al., Nat Methods 2014 | BsmBI |
| `LentiCRISPRv2GFP.gb` | LentiCRISPRv2GFP, Addgene 82416 | Addgene | BsmBI |
| `lentiGuide-Puro.gb` | Addgene 52963 | Sanjana et al., Nat Methods 2014, doi:10.1038/nmeth.3047 | BsmBI |
| `pX330.gb` | pX330-U6-Chimeric_BB-CBh-hSpCas9, Addgene 42230 | Cong et al., Science 2013, doi:10.1126/science.1231143 | BbsI |
| `pU6-BbsI_CBh-Cas9-T2A-mCherry-H1.gb` | Addgene 64217 (U6 for sgRNA, H1 for shRNA) | Chu et al., Nat Biotechnol 2015 | BbsI |
| `phH1-gRNA.gb` | Addgene 53186 (H1 promoter) | Kabadi et al., Nucleic Acids Res 2014, doi:10.1093/nar/gku749 | BbsI |
| `pU6-sgRNA_EF1Alpha-puro-T2A-BFP.gb` | Addgene 60955 | Gilbert et al., Cell 2014 (Weissman lab) | BstXI + BlpI |
| `pMCB306.gb` | Addgene 89360 (mouse U6) | Han et al., Nat Biotechnol 2017, doi:10.1038/nbt.3834 | BstXI + BlpI |

The H1 and U6 reference sequences used to find promoters (`PyShinyApp/promoters.py`) come from these files.
The files for Addgene plasmids 42230, 52963, 53186, 60955 and 64217 were downloaded from Addgene in
November 2026; the first line (LOCUS) of some carries an Addgene sequence number instead of the plasmid name.

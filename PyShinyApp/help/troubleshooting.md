# Troubleshooting

**"BsmBI (site CGTCTC) was found N time(s) in this vector, but exactly 2 are needed"**
The enzyme you chose does not cut your vector exactly twice. The message lists the
supported enzymes that do. If none do, your vector may use a different enzyme or have
extra internal sites; domesticate the extra sites or contact us.

**"The vector GenBank file is not marked circular"**
Only circular plasmids are supported. Export the map as circular from SnapGene or Benchling.

**"The cloning site spans the origin"**
The stuffer crosses position 1 of the map. Move the origin elsewhere in your plasmid editor,
re-export the GenBank file, and try again.

**"Could not confirm the sgRNA orientation"**
The right-hand overhang is not the start of the Cas9 scaffold (`GTTT`) on either strand, so
the app assumed the sgRNA reads forward. Check the overhangs shown under **Cloning site**
and set the orientation under **⚙ Advanced options** if needed.

**"sgRNA ... contains a BsmBI recognition site"**
The enzyme would cut inside this guide during cloning, so it will likely fail. Choose a
different guide.

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


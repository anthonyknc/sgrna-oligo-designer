"""
Find the Pol III promoter that drives the sgRNA in a vector.

The promoter tells us where the sgRNA slot is and which way it is transcribed, so the
cloning site can be located from it instead of only from the enzyme sites. Detection is
tried in this order:

1. Annotation: a feature labelled U6, H1, 7SK or "Pol III promoter" (any spelling).
2. Sequence: if no annotation is found, the 3' end of the promoter is looked for by local
   alignment against reference sequences (human U6, mouse U6, H1).

7SK is recognised from annotations only (no verified reference sequence yet).
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from functools import lru_cache

from Bio import Align
from Bio.Seq import Seq

# 3' ends of the promoters (the ~100 bp that contain the TATA box and PSE, which are the
# best conserved part). Taken from ipUSEPR (human U6), pMCB306 (mouse U6, the part after its
# loxP site) and phH1-gRNA / Addgene 64217 (H1, identical in both); the trailing CACC of the human U6 is left out because vectors differ there.
REFERENCES = {
    "human U6": "AGTTTGCAGTTTTAAAATTATGTTTTAAAATGGACTATCATATGCTTACCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAAC",
    "H1": "GGGACAGGGGAGTGGCGCCCTGCAATATTTGCATGTCGCTATGTGTTCTGGGAAATCACCATAAACGTGAAATGTCTTTGGATTTGGGAATCTTATAAGTTCTGTATGAGACCACT",
    "mouse U6": "CAGCACAAAAGGAAACTCACCCTAACTGTAAAGTAATTGTGTGTTTTGAGACTATAAGTATCCCTTGGAGAACCACCTTGTTG",
}
MIN_IDENTITY = 0.85         # fraction of the reference that must match
MAX_GAP = 150               # a cloning site this far (or less) after the promoter end counts as "downstream"

# annotation labels: (regex, name). Checked in order, on label/note/name/standard_name.
_ANNOTATION_NAMES = [
    (re.compile(r"(?<![a-z0-9])(h|human)\s*[-_ ]?u6(?![a-z0-9])", re.I), "human U6"),
    (re.compile(r"(?<![a-z0-9])(m|mouse|murine)\s*[-_ ]?u6(?![a-z0-9])", re.I), "mouse U6"),
    (re.compile(r"(?<![a-z0-9])u6(?![a-z0-9])", re.I), "U6"),
    (re.compile(r"(?<![a-z0-9])(h|human)?\s*[-_ ]?h1\s*(promoter|prom)(?![a-z0-9])", re.I), "H1"),
    (re.compile(r"(?<![a-z0-9])(h|human)?\s*7sk(?![a-z0-9])", re.I), "7SK"),
    (re.compile(r"pol\s*-?(iii|3)\s*(promoter|prom)", re.I), "Pol III promoter"),
]
_SKIP_TYPES = {"primer_bind", "primer", "CDS", "gene", "source", "mRNA", "misc_binding", "protein_bind"}


@dataclass(frozen=True)
class Pol3Promoter:
    """A promoter on the record. ``start``/``end`` are 0-based, ``end`` exclusive; ``strand`` is
    +1 or -1. ``downstream_end`` is the end of the promoter in the direction of transcription."""
    name: str
    start: int
    end: int
    strand: int
    source: str                 # 'annotation' or 'sequence'
    identity: float = 1.0       # fraction of the reference that matched (sequence only)

    def in_frame(self, n: int, orientation: str) -> int:
        """Position of the promoter's 3' end in the sense frame ('forward' or 'reverse')."""
        forward = self.strand == +1
        if orientation == "forward":
            return self.end if forward else -1          # reverse-strand promoter: not in this frame
        return n - self.start if not forward else -1

    @property
    def label(self) -> str:
        how = "annotated" if self.source == "annotation" else f"found by sequence, {self.identity:.0%} match"
        return f"{self.name} promoter ({how})"


def _annotation_name(feature) -> "str | None":
    if feature.type in _SKIP_TYPES:
        return None
    q = feature.qualifiers
    text = " ".join(" ".join(q.get(k, [])) for k in ("label", "note", "standard_name", "name", "gene"))
    if re.search(r"primer|seq\s*f|seq\s*r|-f$|-r$", text, re.I) and feature.type != "promoter":
        return None
    for pattern, name in _ANNOTATION_NAMES:
        if pattern.search(text):
            return name
    return None


def _from_annotations(record) -> list[Pol3Promoter]:
    found = []
    for f in record.features:
        name = _annotation_name(f)
        if name is None:
            continue
        loc = f.location
        found.append(Pol3Promoter(name, int(loc.start), int(loc.end),
                                  -1 if loc.strand == -1 else +1, "annotation"))
    return found


@lru_cache(maxsize=1)
def _aligner() -> Align.PairwiseAligner:
    a = Align.PairwiseAligner(mode="local")
    a.match_score, a.mismatch_score = 2, -3
    a.open_gap_score, a.extend_gap_score = -6, -3
    return a


def _from_sequence(seq: str) -> list[Pol3Promoter]:
    found, n = [], len(seq)
    for name, ref in REFERENCES.items():
        for strand in (+1, -1):
            target = seq if strand == +1 else str(Seq(seq).reverse_complement())
            aln = _aligner().align(target, ref)[0]
            blocks = aln.aligned[1]                       # aligned parts of the reference
            matched_ref = sum(int(b - a) for a, b in blocks)
            ident = aln.counts().identities / len(ref)
            if matched_ref < MIN_IDENTITY * len(ref) or ident < MIN_IDENTITY:
                continue
            t0, t1 = int(aln.aligned[0][0][0]), int(aln.aligned[0][-1][1])
            start, end = (t0, t1) if strand == +1 else (n - t1, n - t0)
            found.append(Pol3Promoter(name, start, end, strand, "sequence", round(ident, 3)))
    return found


def find_pol3_promoters(record) -> list[Pol3Promoter]:
    """Pol III promoters in a vector, from annotations if there are any, else from the sequence."""
    found = _from_annotations(record)
    if not found:
        found = _from_sequence(str(record.seq).upper())
    return sorted(found, key=lambda p: p.start)

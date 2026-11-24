"""
sgrna_oligo_designer.py
------------------------
Core logic for the sgRNA Oligo Designer web app (see app.py).

The user picks a vector and a Type IIS restriction enzyme. The vector is then
digested *in silico*: the two recognition sites are located, the cut positions
and the 4-nt overhangs are read off the vector, and the annealing oligos and the
final plasmid are derived from those real cut ends. Nothing about the overhangs
or the stuffer has to be typed in.

Worked example (ipUSEPR, BsmBI)::

    ...CACC G GAGACGTGACCGTCTCT GTTT AAGAGC...      vector
       ^^^^                     ^^^^
       left overhang            right overhang       (read from the BsmBI cuts)

    sense     5'-CACC + [G] + sgRNA-3'               overhang + core
    antisense 5'-rc(core + GTTT)-3' = AAAC + rc(sgRNA) + C

This module has no command-line interface; the Shiny app is the only entry point.
"""

import copy
import csv
import re
from collections import Counter
from dataclasses import dataclass, replace as _replace
from functools import lru_cache
from pathlib import Path

from promoters import MAX_GAP as _PROMOTER_GAP, Pol3Promoter, find_pol3_promoters


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# IDT bulk upload template (matches the IDT spec exactly)
IDT_COLUMNS    = ["Name", "Sequence", "Scale", "Purification"]
IDT_SCALE      = "25nm"
IDT_PURIF      = "STD"
IDT_COL_WIDTHS = {"A": 12.54, "B": 41.54, "C": 10.18, "D": 18.82}

# The enzymes listed first in the app, chosen from Addgene's CRISPR plasmids. BsmBI/Esp3I and
# BbsI/BpiI are isoschizomer pairs (same site and cut), so each pair is one entry. Any other
# eligible enzyme can be chosen under "Other"; see all_enzymes(). Cut geometry always comes
# from Bio.Restriction.
ENZYMES = ("BsmBI", "BbsI", "BsaI", "SapI")

# An enzyme is offered only if it is commercially available, recognises a plain ACGT site of
# at least MIN_SITE bases, is non-palindromic, cuts outside its site with a single cut on
# the downstream side, and leaves a 5' overhang of MIN_OVERHANG to MAX_OVERHANG bases
# (1-2 nt overhangs ligate poorly and are not used for sgRNA cloning).
MIN_SITE, MIN_OVERHANG, MAX_OVERHANG = 4, 3, 5

# How the 5' G required by the human U6 promoter is added to the spacer:
#   always    - always prepend G (standard Zhang-lab / Broad convention)
#   if_needed - prepend G only when the spacer does not already start with G
#   never     - never prepend G (non-U6 promoters)
G_MODES = ("always", "if_needed", "never")

# Which way the sgRNA is transcribed relative to the GenBank file.
ORIENTATIONS = ("auto", "forward", "reverse")

# The first four bases of the SpCas9 scaffold. The right-hand overhang equals
# this in essentially every SpCas9 vector, which lets us detect orientation.
SCAFFOLD_START = "GTTT"

# Start of the Cas9 scaffold, allowing the common variants (GTTTTAGAGC, improved GTTTAAGAGC).
SCAFFOLD_RE = re.compile("GTTT[ACGT]AGAGC")
SPACER_LEN  = 20      # assumed length of the vector's placeholder spacer when suggesting flanks

# Largest sgRNA table accepted in one run (keeps the web app responsive).
MAX_SGRNAS = 1000

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


_IUPAC_RE = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]", "S": "[CG]",
             "W": "[AT]", "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]", "H": "[ACT]",
             "V": "[ACG]", "N": "[ACGT]"}
_IUPAC_COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "R": "Y", "Y": "R", "S": "S", "W": "W",
               "K": "M", "M": "K", "B": "V", "V": "B", "D": "H", "H": "D", "N": "N"}


def _site_regexes(site: str) -> tuple[str, str]:
    """Regexes matching an (IUPAC) recognition site on the top and on the bottom strand."""
    site = site.upper()
    return (
        "".join(_IUPAC_RE[c] for c in site),
        "".join(_IUPAC_RE[_IUPAC_COMP[c]] for c in reversed(site)),
    )


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence (ACGT only)."""
    return seq.translate(_COMPLEMENT)[::-1]


# ---------------------------------------------------------------------------
# Restriction enzymes and in-silico digestion
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Enzyme:
    """Cut geometry of a Type IIS enzyme that leaves a 5' overhang."""
    name: str
    site: str           # recognition sequence, top strand, ACGT only
    cut_offset: int     # top-strand cut, measured from the start of the site
    overhang_len: int   # length of the 5' overhang


def _from_biopython(enz) -> "Enzyme | None":
    """Return the cut geometry of a Bio.Restriction enzyme, or None if it is not eligible."""
    site = str(enz.site).upper()
    if (not enz.suppl or set(site) - set("ACGT") or len(site) < MIN_SITE
            or enz.is_palindromic() or not enz.is_5overhang()
            or enz.fst5 is None or enz.fst5 < len(site)
            or enz.scd5 is not None or enz.scd3 is not None
            or not MIN_OVERHANG <= -enz.ovhg <= MAX_OVERHANG):
        return None
    return Enzyme(name=enz.__name__, site=site, cut_offset=enz.fst5, overhang_len=-enz.ovhg)


@lru_cache(maxsize=1)
def all_enzymes() -> dict[str, Enzyme]:
    """Every eligible enzyme in Bio.Restriction's database, keyed by name and sorted."""
    from Bio import Restriction

    found = (_from_biopython(e) for e in Restriction.AllEnzymes)
    return {e.name: e for e in sorted((e for e in found if e), key=lambda e: e.name.lower())}


def get_enzyme(name: str) -> Enzyme:
    """Look up an eligible enzyme by name (see :func:`all_enzymes`)."""
    enzymes = all_enzymes()
    if name not in enzymes:
        raise ValueError(
            f"Unsupported enzyme '{name}'. It must be a commercially available Type IIS enzyme with a "
            f"plain recognition site that leaves a {MIN_OVERHANG}-{MAX_OVERHANG} nt 5' overhang."
        )
    return enzymes[name]


def _find_sites(seq: str, enzyme: Enzyme) -> list[tuple[int, int]]:
    """
    Find every recognition site on a circular sequence, on both strands.

    Returns sorted ``(start, strand)`` tuples. ``start`` is the 0-based index of
    the first base of the match on the top strand; strand is +1 when the site
    reads 5'->3' on the top strand and -1 when its reverse complement does.
    """
    n, site_len = len(seq), len(enzyme.site)
    extended    = seq + seq[:site_len - 1]          # catch sites spanning the origin
    hits = []
    for strand, pattern in ((+1, enzyme.site), (-1, reverse_complement(enzyme.site))):
        pos = 0
        while (i := extended.find(pattern, pos)) != -1 and i < n:
            hits.append((i, strand))
            pos = i + 1
    return sorted(hits)


def _overhang_start(start: int, strand: int, enzyme: Enzyme) -> int:
    """Top-strand index where the overhang left by one site begins."""
    if strand == +1:
        return start + enzyme.cut_offset
    return start + len(enzyme.site) - enzyme.cut_offset - enzyme.overhang_len


@dataclass(frozen=True)
class CloningSite:
    """
    Result of digesting a vector in silico. All coordinates are 0-based and refer
    to the sequence in the *sense frame*, i.e. the vector read in the direction
    the sgRNA is transcribed (the GenBank sequence itself if ``orientation`` is
    'forward', its reverse complement if 'reverse').
    """
    enzyme:         str
    enzyme_site:    str
    orientation:    str     # 'forward' or 'reverse'
    left_overhang:  str     # sense-frame sequence, 5' side of the insert
    right_overhang: str     # sense-frame sequence, 3' side of the insert
    replace_start:  int     # the bases [replace_start, replace_end) are replaced
    replace_end:    int
    left_context:   str     # backbone bases just outside each overhang, used to
    right_context:  str     # detect enzyme sites created at the junctions
    vector_length:  int
    warnings:       tuple = ()
    # --- two-enzyme / general geometry (defaults describe one 5'-overhang enzyme on both sides)
    enzyme_b:         str = ""      # second enzyme, if the two ends are cut by different enzymes
    site_patterns:    tuple = ()    # regexes for every recognition site; empty = derive from enzyme_site
    left_extra:       str = ""      # constant bases between the left overhang and the spacer
    right_extra:      str = ""      # constant bases between the spacer and the right overhang
    left_top_trim:    int = 0       # bases of the assembled insert missing from each oligo's ends;
    left_bottom_trim: "int | None" = None    # None = the overhang length (a 5' overhang end)
    right_top_trim:   "int | None" = None
    right_bottom_trim: int = 0
    promoter:         str = ""      # the Pol III promoter just upstream of the site, e.g. 'human U6 promoter (annotated)'

    @property
    def removed_length(self) -> int:
        """Bases between the two overhangs that are replaced by the insert."""
        return self.replace_end - self.replace_start

    @property
    def left_kind(self) -> str:
        """Overhang type at the left end of the insert: 5' or 3'."""
        return "3'" if self.left_top_trim > 0 else "5'"

    @property
    def right_kind(self) -> str:
        """Overhang type at the right end of the insert: 5' or 3'."""
        return "3'" if self.right_bottom_trim > 0 else "5'"

    @property
    def enzymes(self) -> str:
        """Display name of the enzyme(s), e.g. 'BsmBI' or 'BstXI + BlpI'."""
        return f"{self.enzyme} + {self.enzyme_b}" if self.enzyme_b else self.enzyme


def _crosses_origin_message(label: str, noun: str, removed: int, other: int) -> str:
    return (
        f"The two {label} {noun} split the plasmid into a {removed} bp piece between the "
        f"{noun} and a {other} bp piece that crosses the origin (position 1) of the map. "
        "The shorter piece is taken to be the stuffer, and a stuffer that crosses the origin is "
        "not supported. Move the origin elsewhere in your plasmid editor, re-export the "
        "GenBank file, and try again. If that does not fit your vector, check that "
        f"{label} is what is used to clone into it (it may just happen to cut there by chance), "
        "or that the stuffer is not the longer piece, which is not supported yet."
    )


def _locate(seq: str, enzyme: Enzyme) -> tuple[int, int]:
    """
    Digest ``seq`` in silico and return the start of the left and right overhang.

    The vector must have exactly two sites. The piece removed is the shorter of
    the two circular fragments (the stuffer); it must not span the origin.
    Both layouts work: recognition sites inside the removed piece (lost on
    cutting, e.g. lentiCRISPR) or left behind in the backbone.
    """
    n, k  = len(seq), enzyme.overhang_len
    hits  = _find_sites(seq, enzyme)
    if len(hits) != 2:
        others = [name for name, e in all_enzymes().items()
                  if name != enzyme.name and len(_find_sites(seq, e)) == 2]
        # familiar enzymes first, then longest (most specific) site first: a shorter site that is
        # contained in a longer one (e.g. EarI inside SapI) cuts there too but is rarely the intended enzyme
        others.sort(key=lambda n: (n not in ENZYMES, -len(all_enzymes()[n].site), n.lower()))
        shown = ", ".join(others[:10]) + (f" and {len(others) - 10} more" if len(others) > 10 else "")
        hint = (f" This vector has exactly two sites for: {shown}." if others
                else " None of the other supported enzymes cuts it exactly twice.")
        raise ValueError(
            f"{enzyme.name} (site {enzyme.site}) was found {len(hits)} time(s) in this vector, "
            f"but exactly 2 are needed, one on each side of the stuffer.{hint}"
        )

    x0, y0 = sorted(_overhang_start(i, s, enzyme) for i, s in hits)
    origin_msg = (
        f"The {enzyme.name} cloning site spans the origin (position 1) of the plasmid map. "
        "Move the origin elsewhere in your plasmid editor, re-export the GenBank file, and try again."
    )
    if x0 < 0 or y0 + k > n:
        raise ValueError(origin_msg)
    if y0 < x0 + k:
        raise ValueError(
            f"The two {enzyme.name} cuts overlap, so no insert can be placed between them."
        )
    removed, other_piece = y0 - (x0 + k), n - (y0 + k - x0)
    if removed > other_piece:
        raise ValueError(_crosses_origin_message(enzyme.name, "sites", removed, other_piece))
    return x0, y0


def _site_from(seq: str, enzyme: Enzyme, x0: int, y0: int, orientation: str,
               warnings: tuple = ()) -> CloningSite:
    k, ctx = enzyme.overhang_len, len(enzyme.site) - 1
    left, right = seq[x0:x0 + k], seq[y0:y0 + k]
    notes = list(warnings)
    if left == right:
        notes.append(
            f"Both overhangs are identical ({left}), so the insert can ligate in either "
            "orientation. Check the vector."
        )
    elif left == reverse_complement(left) or right == reverse_complement(right):
        notes.append(
            f"An overhang ({left if left == reverse_complement(left) else right}) is palindromic, "
            "so vector ends can ligate to each other and empty vectors are likely."
        )
    return CloningSite(
        enzyme=enzyme.name,
        enzyme_site=enzyme.site,
        orientation=orientation,
        left_overhang=seq[x0:x0 + k],
        right_overhang=seq[y0:y0 + k],
        replace_start=x0 + k,
        replace_end=y0,
        left_context=seq[max(0, x0 - ctx):x0],
        right_context=seq[y0 + k:y0 + k + ctx],
        vector_length=len(seq),
        warnings=tuple(notes),
        site_patterns=_site_regexes(enzyme.site),
    )


# ---------------------------------------------------------------------------
# Two-enzyme cloning (the two ends are cut by different enzymes)
# ---------------------------------------------------------------------------

@lru_cache(maxsize=1)
def two_enzyme_candidates() -> dict:
    """
    Every commercially available enzyme (from Bio.Restriction) that can serve as one end of a
    two-enzyme cloning site: it leaves a 5' or 3' overhang and cuts in one place. Palindromic
    and degenerate sites (e.g. EcoRI, AgeI, BstXI, BlpI) are fine here. Keyed by name.
    """
    from Bio import Restriction

    found = {
        e.__name__: e for e in Restriction.AllEnzymes
        if e.suppl and e.ovhg and e.fst5 is not None and e.scd5 is None and e.scd3 is None
    }
    return dict(sorted(found.items(), key=lambda kv: kv[0].lower()))


def _cuts(seq: str, enz) -> list[tuple[int, int]]:
    """
    Cuts of ``enz`` on a circular sequence as (top, bottom) boundary positions: the top strand
    is cut between bases top-1 and top, the bottom strand between bases bottom-1 and bottom
    (both in top-strand coordinates). A 5' overhang has bottom > top, a 3' overhang bottom < top.
    """
    from Bio.Seq import Seq

    n = len(seq)
    out = []
    for p in enz.search(Seq(seq), linear=False):
        top = (p - 1) % n
        out.append((top, top - enz.ovhg))
    return out


def single_cutters(vector_record) -> list[str]:
    """Names of the enzymes usable in two-enzyme mode that cut this vector exactly once."""
    seq = str(vector_record.seq).upper()
    return [name for name, enz in two_enzyme_candidates().items() if len(_cuts(seq, enz)) == 1]


def _locate_two(seq: str, enz_a, enz_b):
    """
    Digest ``seq`` in silico with two different enzymes that each cut once. Returns the left
    and right cut as (enzyme, top, bottom, lowest, highest) tuples, where lowest/highest bound
    the overhang. The piece between the cuts is the stuffer: the shorter of the two circular
    pieces, which must not cross the origin.
    """
    n, cuts = len(seq), []
    for enz in (enz_a, enz_b):
        found = _cuts(seq, enz)
        if len(found) != 1:
            raise ValueError(
                f"{enz.__name__} cuts this vector {len(found)} time(s), but each of the two enzymes "
                "must cut it exactly once."
            )
        top, bottom = found[0]
        cuts.append((enz, top, bottom, min(top, bottom), max(top, bottom)))
    left, right = sorted(cuts, key=lambda c: c[3])
    if left[3] < 0 or right[4] > n:
        raise ValueError(
            f"The {left[0].__name__}/{right[0].__name__} cloning site spans the origin (position 1) "
            "of the plasmid map. Move the origin elsewhere in your plasmid editor, re-export the "
            "GenBank file, and try again."
        )
    if left[4] > right[3]:
        raise ValueError(
            f"The {left[0].__name__} and {right[0].__name__} cuts overlap, so no insert can be "
            "placed between them."
        )
    removed, other = right[3] - left[4], n - (right[4] - left[3])
    if removed > other:
        raise ValueError(_crosses_origin_message(
            f"{left[0].__name__}/{right[0].__name__}", "cuts", removed, other))
    return left, right


def _suggest_flanks(seq: str, left_max: int, right_min: int, right_max: int):
    """
    Guess the constant bases around the spacer from the vector itself. Returns
    ``(left_extra, right_extra, scaffold_found)``.

    The Cas9 scaffold start (GTTT.AGAGC) marks the 3' end of the spacer slot, so the bases from
    its first G up to the right overhang are the right flank. The slot is assumed to be
    ``SPACER_LEN`` bases long, so the bases before it (after the left overhang) are the left flank.
    """
    m = SCAFFOLD_RE.search(seq, left_max, min(len(seq), right_max + 10))
    if m is None or m.start() >= right_max:
        return "", "", False
    right = seq[m.start():right_min] if m.start() < right_min else ""
    left = seq[left_max:m.start() - SPACER_LEN] if m.start() - SPACER_LEN >= left_max else ""
    return left, right, True


def _clean_flank(value: str, label: str) -> str:
    flank = "".join(value.split()).upper()
    if set(flank) - set("ACGT"):
        raise ValueError(f"{label} may contain only A, C, G and T.")
    return flank


def _site_two(seq, left, right, orientation, left_extra, right_extra):
    """Build a CloningSite from a two-enzyme digest. Returns (site, scaffold_found)."""
    (ea, ta, ba, ma, Ma), (eb, tb, bb, mb, Mb) = left, right
    left_ov, right_ov = seq[ma:Ma], seq[mb:Mb]
    sug_l, sug_r, found = _suggest_flanks(seq, Ma, mb, Mb)
    l_extra = sug_l if left_extra is None else _clean_flank(left_extra, "The 5' flanking bases")
    r_extra = sug_r if right_extra is None else _clean_flank(right_extra, "The 3' flanking bases")

    notes = []
    if left_ov == right_ov:
        notes.append(f"Both overhangs are identical ({left_ov}), so the insert can ligate in either "
                     "orientation. Check the vector.")
    pal = [ov for ov in (left_ov, right_ov) if ov == reverse_complement(ov)]
    if pal and left_ov != right_ov:
        notes.append(f"An overhang ({pal[0]}) is palindromic, so vector ends can ligate to each other "
                     "and empty vectors are likely.")
    if left_extra is None and right_extra is None and not found:
        notes.append("The Cas9 scaffold start was not found between the cuts, so the constant bases "
                     "around the spacer could not be worked out. Enter them under Advanced options "
                     "(as in your vector's oligo design instructions).")

    pats = tuple(dict.fromkeys(
        pat for e in (ea, eb) for pat in _site_regexes(str(e.site))))
    ctx = max(len(str(ea.site)), len(str(eb.site))) - 1
    site = CloningSite(
        enzyme=ea.__name__,
        enzyme_b=eb.__name__,
        enzyme_site=f"{ea.site} / {eb.site}",
        orientation=orientation,
        left_overhang=left_ov,
        right_overhang=right_ov,
        replace_start=Ma,
        replace_end=mb,
        left_context=seq[max(0, ma - ctx):ma],
        right_context=seq[Mb:Mb + ctx],
        vector_length=len(seq),
        warnings=tuple(notes),
        site_patterns=pats,
        left_extra=l_extra,
        right_extra=r_extra,
        left_top_trim=ta - ma,
        left_bottom_trim=ba - ma,
        right_top_trim=Mb - tb,
        right_bottom_trim=Mb - bb,
    )
    return site, found


def _find_two_enzyme_site(seq, rc_seq, name_a, name_b, orientation, left_extra, right_extra,
                          promoters=()):
    candidates = two_enzyme_candidates()
    for name in (name_a, name_b):
        if name not in candidates:
            raise ValueError(
                f"Unsupported enzyme '{name}'. It must be commercially available and leave a 5' or "
                "3' overhang."
            )
    if name_a == name_b:
        raise ValueError("Choose two different enzymes (or use one-enzyme mode).")
    ea, eb = candidates[name_a], candidates[name_b]

    def build(frame, tag):
        left, right = _locate_two(frame, ea, eb)
        return _site_two(frame, left, right, tag, left_extra, right_extra)

    if orientation == "forward":
        return build(seq, "forward")[0]
    if orientation == "reverse":
        return build(rc_seq, "reverse")[0]
    fwd, fwd_found = build(seq, "forward")
    if fwd_found:
        return fwd
    rev, rev_found = build(rc_seq, "reverse")
    if rev_found:
        return rev
    n = len(seq)
    pf, pr = (_promoter_before(x, list(promoters), n) is not None for x in (fwd, rev))
    if pf != pr:                                      # no scaffold, but the promoter shows the way
        return fwd if pf else rev
    return CloningSite(**{**fwd.__dict__, "warnings": fwd.warnings + (
        "Could not confirm the sgRNA orientation: the Cas9 scaffold start was not found on either "
        "strand. Assuming forward; set the orientation under Advanced options if this is wrong.",)})


def _promoter_before(site: CloningSite, promoters: list, n: int) -> "Pol3Promoter | None":
    """The promoter that ends just upstream of the site in its own sense frame, if any."""
    for p in promoters:
        end = p.in_frame(n, site.orientation)
        if end >= 0 and -8 <= site.replace_start - end <= _PROMOTER_GAP:
            return p
    return None


def _annotate(site: CloningSite, promoters: list, n: int) -> CloningSite:
    p = _promoter_before(site, promoters, n)
    return _replace(site, promoter=p.label) if p else site


def _looks_like_scaffold_start(right_overhang: str) -> bool:
    """
    True if a right-hand overhang is the start of the Cas9 scaffold (GTTT...). Needs at
    least 3 bases to be informative; longer overhangs are compared on their first 4 bases.
    """
    n = min(len(right_overhang), len(SCAFFOLD_START))
    return len(right_overhang) >= 3 and right_overhang[:n] == SCAFFOLD_START[:n]


def find_cloning_site(
    vector_record,
    enzyme_name: str,
    orientation: str = "auto",
    *,
    enzyme_b: "str | None" = None,
    left_extra: "str | None" = None,
    right_extra: "str | None" = None,
) -> CloningSite:
    """
    Digest the vector in silico with ``enzyme_name`` and describe the cloning site.

    Parameters
    ----------
    vector_record : Bio.SeqRecord.SeqRecord
        Circular vector.
    enzyme_name : str
        Any name from :func:`all_enzymes` (the common ones are in :data:`ENZYMES`).
    orientation : str
        'forward', 'reverse', or 'auto'. Auto-detection looks for the SpCas9
        scaffold start (GTTT) in the right-hand overhang on either strand and
        falls back to 'forward' with a warning if neither matches.
    enzyme_b : str, optional
        Two-enzyme mode: a second, different enzyme. ``enzyme_name`` and ``enzyme_b``
        must each cut the vector exactly once (see :func:`single_cutters`), and may be
        any commercially available enzyme (palindromic, degenerate, 3' or 5' overhang).
    left_extra, right_extra : str, optional
        Two-enzyme mode only: the constant bases that sit between the overhangs and the
        spacer (for BstXI/BlpI vectors, e.g. ``TTG`` and ``GTTTAAGAGC``). ``None`` means
        "work them out from the vector" (see :func:`_suggest_flanks`); ``""`` means none.

    Raises
    ------
    ValueError
        If the vector is not circular, does not have exactly two sites for the
        enzyme, or the cloning site spans the plasmid origin.
    """
    if orientation not in ORIENTATIONS:
        raise ValueError(f"orientation must be one of {ORIENTATIONS}, got '{orientation}'.")
    enzyme = None if enzyme_b else get_enzyme(enzyme_name)
    if vector_record.annotations.get("topology") != "circular":
        raise ValueError(
            "The vector GenBank file is not marked circular. Only circular plasmids are supported."
        )

    seq    = str(vector_record.seq).upper()
    rc_seq = reverse_complement(seq)
    n = len(seq)
    promoters = find_pol3_promoters(vector_record)
    if enzyme_b:
        return _annotate(_find_two_enzyme_site(seq, rc_seq, enzyme_name, enzyme_b, orientation,
                                                    left_extra, right_extra, promoters), promoters, n)
    if left_extra is not None or right_extra is not None:
        raise ValueError("Flanking bases can only be given in two-enzyme mode.")
    fx0, fy0 = _locate(seq, enzyme)                  # raises with a helpful message

    def forward(warnings=()):
        return _site_from(seq, enzyme, fx0, fy0, "forward", warnings)

    def reverse():
        rx0, ry0 = _locate(rc_seq, enzyme)
        return _site_from(rc_seq, enzyme, rx0, ry0, "reverse")

    if orientation == "forward":
        return _annotate(forward(), promoters, n)
    if orientation == "reverse":
        return _annotate(reverse(), promoters, n)

    fwd, rev = forward(), reverse()
    scaf_f, scaf_r = (_looks_like_scaffold_start(x.right_overhang) for x in (fwd, rev))
    prom_f, prom_r = (_promoter_before(x, promoters, n) is not None for x in (fwd, rev))
    if scaf_f != scaf_r:
        chosen = fwd if scaf_f else rev
        if prom_f != prom_r and (prom_f != scaf_f):
            chosen = _replace(chosen, warnings=chosen.warnings + (
                "The Cas9 scaffold start and the Pol III promoter point to different orientations. "
                "Check the orientation under Advanced options.",))
        return _annotate(chosen, promoters, n)
    if prom_f != prom_r:                              # the promoter decides when the scaffold cannot
        return _annotate(fwd if prom_f else rev, promoters, n)
    if scaf_f:                                        # both look like a scaffold start
        return _annotate(fwd, promoters, n)
    return _annotate(forward((
        f"Could not confirm the sgRNA orientation: expected the scaffold start ({SCAFFOLD_START}) "
        f"at the start of the right-hand overhang, found {fwd.right_overhang} (forward) / {rev.right_overhang} "
        "(reverse). Assuming forward; set the orientation under Advanced options if this is wrong.",
    )), promoters, n)


def suggest_cloning_sites(vector_record, limit: int = 12) -> list[CloningSite]:
    """
    Cloning sites that sit just downstream of the vector's Pol III promoter, best first.

    One-enzyme sites (enzymes that cut the vector exactly twice) come first, familiar enzymes
    before the rest, then two-enzyme sites (two enzymes that each cut once, the left one right
    after the promoter). Returns an empty list if no promoter is found or nothing fits.
    """
    promoters = find_pol3_promoters(vector_record)
    if not promoters:
        return []
    out: list[CloningSite] = []
    for name in all_enzymes():
        try:
            site = find_cloning_site(vector_record, name)
        except ValueError:
            continue
        if site.promoter:
            out.append(site)
    out.sort(key=lambda x: (x.enzyme not in ENZYMES, x.enzyme.lower()))

    seq = str(vector_record.seq).upper()
    n = len(seq)
    cuts = {name: _cuts(seq, enz) for name, enz in two_enzyme_candidates().items()}
    singles = {name: c[0] for name, c in cuts.items() if len(c) == 1}
    for p in promoters:
        orient = "forward" if p.strand == +1 else "reverse"
        end = p.in_frame(n, orient)
        frame = seq if orient == "forward" else reverse_complement(seq)
        scaffold = SCAFFOLD_RE.search(frame, end, end + 400)
        if scaffold is None:
            continue
        pos = {name: (min(c) if orient == "forward" else n - max(c)) for name, c in singles.items()}
        left = [a for a, pa in pos.items()
                if -8 <= pa - end <= _PROMOTER_GAP and pa <= scaffold.start() - SPACER_LEN + 5]
        right = [b for b, pb in pos.items() if scaffold.start() - 2 <= pb <= scaffold.start() + 40]
        for a in left:
            for b in right:
                try:
                    site = find_cloning_site(vector_record, a, orient, enzyme_b=b)
                except ValueError:
                    continue
                if site.promoter:
                    out.append(site)
    # isoschizomers cut identically: keep one entry per distinct cut
    unique, seen = [], set()
    for x in out:
        key = (x.orientation, x.left_overhang, x.right_overhang, x.replace_start, x.replace_end,
               x.left_top_trim, x.right_bottom_trim)
        if key not in seen:
            seen.add(key)
            unique.append(x)
    out = unique
    return out[:limit]


# ---------------------------------------------------------------------------
# Oligo design
# ---------------------------------------------------------------------------

def _validate_sgrna_seq(seq: str, sgrna_id: str) -> tuple[str, list[str]]:
    """
    Clean and validate a raw sgRNA sequence.

    Removes all whitespace (including copy-paste artifacts) and uppercases.
    Raises ValueError for empty or non-ACGT sequences. Returns the cleaned
    sequence and a list of warnings for a length other than 20 nt and for a
    likely PAM (trailing NGG on a sequence longer than 20 nt).
    """
    seq = "".join(seq.split()).upper()
    warnings: list[str] = []

    if not seq:
        raise ValueError(
            f"sgRNA '{sgrna_id}': sequence is empty. Check that the sequence column is correct."
        )
    invalid = set(seq) - set("ACGT")
    if invalid:
        raise ValueError(
            f"sgRNA '{sgrna_id}': contains non-ACGT characters: {sorted(invalid)}. "
            "Only standard ACGT bases are supported. Remove ambiguity codes, "
            "spaces, or non-DNA characters."
        )
    if len(seq) != 20:
        warnings.append(f"is {len(seq)} nt (expected 20 nt for SpCas9); proceeding anyway.")
    if len(seq) > 20 and seq.endswith("GG"):
        warnings.append(
            f"is {len(seq)} nt and ends in '{seq[-3:]}', which looks like an NGG PAM "
            "included by mistake. Trim to 20 nt if so."
        )
    return seq, warnings


def build_insert_core(sgrna_seq: str, g_mode: str = "always") -> str:
    """Return the spacer with the optional 5' G for U6 transcription added."""
    if g_mode not in G_MODES:
        raise ValueError(f"g_mode must be one of {G_MODES}, got '{g_mode}'.")
    if g_mode == "always" or (g_mode == "if_needed" and not sgrna_seq.startswith("G")):
        return "G" + sgrna_seq
    return sgrna_seq


def design_oligos(core: str, site: CloningSite) -> tuple[str, str]:
    """
    Build the top (sense) and bottom (antisense) oligos for one insert.

    The assembled insert is ``left overhang + left flank + core + right flank + right overhang``
    (the overhangs come from the vector). Each oligo is that strand of the insert, minus the
    bases the vector end already supplies on that strand (recorded in the site's ``*_trim``
    fields). For one 5'-overhang enzyme this gives the classic
    ``sense = left overhang + core`` and ``antisense = rc(core + right overhang)``.
    """
    ins = site.left_overhang + site.left_extra + core + site.right_extra + site.right_overhang
    l_bot = len(site.left_overhang) if site.left_bottom_trim is None else site.left_bottom_trim
    r_top = len(site.right_overhang) if site.right_top_trim is None else site.right_top_trim
    top = ins[site.left_top_trim:len(ins) - r_top]
    bottom = reverse_complement(ins[l_bot:len(ins) - site.right_bottom_trim])
    return top, bottom


def add_cloning_oligos(
    rows: list[dict],
    site: CloningSite,
    seq_col: str = "sgRNA_seq",
    id_col: str  = "sgRNA_id",
    g_mode: str  = "always",
) -> list[dict]:
    """
    Validate each sgRNA and add its cloning oligos, in place.

    Adds to each row: ``insert_core``, ``oligo_sense``, ``oligo_antisense`` and
    ``warnings`` (a list of human-readable strings, empty if all is well). The
    cleaned sequence is written back to ``row[seq_col]``.

    A warning is raised if the spacer (or a junction it creates) contains the
    enzyme's recognition site, because the enzyme would cut the insert during
    cloning.
    """
    if rows and seq_col not in rows[0]:
        raise ValueError(
            f"Column '{seq_col}' not found. Available: {list(rows[0].keys())}"
        )
    patterns = site.site_patterns or _site_regexes(site.enzyme_site)

    for i, row in enumerate(rows):
        sgrna_id = row.get(id_col, f"row {i + 1}")
        seq, notes = _validate_sgrna_seq(row[seq_col], sgrna_id)
        core = build_insert_core(seq, g_mode)

        if "TTTT" in seq:
            notes.append(
                "contains TTTT, which ends transcription from U6 and other RNA polymerase III "
                "promoters, so the sgRNA may be truncated. Choose another guide if your vector "
                "uses such a promoter."
            )
        # Only sites that really lie in the spacer count. The vector's own flanks may legitimately
        # regenerate a site (e.g. BstXI's final G is the guide's first G), so a site that overlaps the
        # spacer by a single base is determined by the vector, not by the guide, and is ignored.
        before = site.left_context + site.left_overhang + site.left_extra
        window = before + core + site.right_extra + site.right_overhang + site.right_context
        s0 = len(before) + (len(core) - len(seq))
        s1 = len(before) + len(core)
        if any(
            min(m.start() + len(m.group(1)), s1) - max(m.start(), s0) >= 2
            for pat in patterns for m in re.finditer(f"(?=({pat}))", window)
        ):
            notes.append(
                f"contains a {site.enzymes} recognition site ({site.enzyme_site}); the enzyme "
                "would cut this insert during cloning, so this guide will likely fail."
            )

        row[seq_col]           = seq
        row["insert_core"]     = core
        row["oligo_sense"], row["oligo_antisense"] = design_oligos(core, site)
        row["warnings"]        = [f"sgRNA '{sgrna_id}' {n}" for n in notes]
    return rows


_NAME_OK  = re.compile(r"[A-Za-z0-9_][A-Za-z0-9_.\-]*")      # no leading - or . (spreadsheet formulas, hidden files)
_NAME_BAD = re.compile(r"[^A-Za-z0-9_.\-]+")


def sanitize_name(name: str, default: str = "unnamed") -> str:
    """
    Make a string safe to use in a file name or a GenBank LOCUS name: runs of any
    character other than letters, digits, '_', '.' and '-' become '_', and leading or
    trailing '.' / '_' are dropped (so '../x' and '/abs/x' can never leave a directory).
    """
    return _NAME_BAD.sub("_", str(name).strip()).strip("._") or default


def validate_ids(rows: list[dict], id_col: str = "sgRNA_id") -> None:
    """
    Check the sgRNA ID column and strip stray whitespace around each ID (in place).

    IDs end up in file names, GenBank LOCUS names and the IDT order form, so they may
    contain only letters, digits, '_', '-' and '.'. In particular this rules out spaces
    (which break GenBank LOCUS lines) and a leading '=' (which Excel runs as a formula).
    Raises ValueError for a missing column, invalid characters, or duplicate IDs.
    """
    if rows and id_col not in rows[0]:
        raise ValueError(
            f"Column '{id_col}' not found. Available: {list(rows[0].keys())}"
        )
    for r in rows:
        r[id_col] = str(r[id_col]).strip()
    bad = [r[id_col] for r in rows if not _NAME_OK.fullmatch(r[id_col])]
    if bad:
        raise ValueError(
            f"Invalid sgRNA ID(s): {bad[:5]}{' ...' if len(bad) > 5 else ''}. IDs may contain only "
            "letters, digits, underscores (_), hyphens (-) and dots (.), must start with a letter, digit "
            "or underscore, and cannot contain spaces. Rename them in your table, for example by "
            "replacing spaces with underscores."
        )
    # Case-insensitive: SG1 and sg1 would give files that overwrite each other on macOS and Windows.
    counts = Counter(r[id_col].lower() for r in rows)
    duplicates = sorted({r[id_col] for r in rows if counts[r[id_col].lower()] > 1})
    if duplicates:
        raise ValueError(
            f"Duplicate sgRNA IDs found: {duplicates} (IDs are compared ignoring upper/lower case). "
            "Each sgRNA must have a unique identifier."
        )


def build_idt_rows(rows: list[dict], id_col: str = "sgRNA_id") -> list[dict]:
    """
    Reshape oligo rows into IDT bulk upload format. Each sgRNA produces two
    interleaved rows:
        <sgRNA_id>_s | oligo_sense     | 25nm | STD
        <sgRNA_id>_a | oligo_antisense | 25nm | STD
    """
    idt_rows = []
    for r in rows:
        for suffix, key in (("s", "oligo_sense"), ("a", "oligo_antisense")):
            idt_rows.append({
                "Name":         f"{r[id_col]}_{suffix}",
                "Sequence":     r[key],
                "Scale":        IDT_SCALE,
                "Purification": IDT_PURIF,
            })
    return idt_rows


# ---------------------------------------------------------------------------
# GenBank plasmid construction
# ---------------------------------------------------------------------------

def read_vector_gb(path: Path):
    """Parse a vector GenBank file into a Bio.SeqRecord.SeqRecord."""
    from Bio import SeqIO
    try:
        return SeqIO.read(path, "genbank")
    except Exception as e:
        msg = str(e)
        if "More than one record" in msg:
            raise ValueError(
                "This GenBank file contains more than one record. Upload a file with a single plasmid."
            ) from e
        raise ValueError(
            "Could not read this file as a GenBank file. Export the map from SnapGene or "
            "Benchling in GenBank (.gb) format."
        ) from e


def _remap_feature(feature, r0: int, r1: int, new_len: int):
    """
    Return a copy of ``feature`` remapped for replacing bases [r0, r1) with
    ``new_len`` new bases, or None if the feature disappears.

    - Entirely before the region: unchanged. Entirely after: shifted.
    - Spanning the region: end shifted.
    - Partly overlapping: clipped to the part that survives.
    - Entirely inside: removed.
    """
    from Bio.SeqFeature import CompoundLocation, FeatureLocation

    loc = feature.location
    delta, parts, changed = new_len - (r1 - r0), [], False
    for p in loc.parts:
        s, e = int(p.start), int(p.end)
        if e <= r0:
            parts.append(FeatureLocation(s, e, p.strand))
        elif s >= r1:
            parts.append(FeatureLocation(s + delta, e + delta, p.strand))
        else:
            changed = True
            ns = s if s < r0 else r0 + new_len
            ne = e + delta if e > r1 else r0
            if ne > ns:
                parts.append(FeatureLocation(ns, ne, p.strand))
    if not parts:
        return None

    new = copy.deepcopy(feature)
    if changed:
        new.qualifiers.pop("translation", None)    # no longer matches the edited sequence
    new.location = parts[0] if len(parts) == 1 else CompoundLocation(parts, operator=loc.operator)
    return new


_RC_KEEP = dict(id=True, name=True, description=True, features=True,
                annotations=True, letter_annotations=True, dbxrefs=True)


def _common_prefix_len(a: str, b: str) -> int:
    n = 0
    for x, y in zip(a, b):
        if x != y:
            break
        n += 1
    return n


def insert_sgrna(vector_record, site: CloningSite, core: str, sgrna_seq: str, sgrna_id: str):
    """
    Replace the region between the two overhangs with ``core`` and annotate the sgRNA.

    Features inside the replaced region are removed, features that overlap it
    are clipped, and everything downstream is shifted. A new ``misc_feature``
    labelled with ``sgrna_id`` marks the spacer (without the added 5' G).

    Parameters
    ----------
    vector_record : Bio.SeqRecord.SeqRecord
        The vector the ``site`` was computed from (not modified).
    site : CloningSite
        From :func:`find_cloning_site`.
    core : str
        Spacer with the optional 5' G, from :func:`build_insert_core`.
    sgrna_seq : str
        The cleaned 20-mer; must be the end of ``core``.
    sgrna_id : str
        Used for the feature label and the record ID.
    """
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature

    if not core.endswith(sgrna_seq):
        raise ValueError("core must end with the sgRNA sequence.")
    reverse = site.orientation == "reverse"
    frame   = vector_record.reverse_complement(**_RC_KEEP) if reverse else vector_record

    seq = str(frame.seq)
    if len(seq) != site.vector_length:
        raise ValueError("The cloning site was computed from a different vector.")
    r0, r1 = site.replace_start, site.replace_end
    g_len  = len(core) - len(sgrna_seq)
    prefix = site.left_extra + core[:g_len]          # constant bases before the spacer
    suffix = site.right_extra                         # constant bases after the spacer
    new_part = prefix + sgrna_seq + suffix
    # Bases that are already identical at the edges of the replaced region (the vector's own
    # 5' G, flanks, ...) stay untouched, so their annotations survive.
    old = seq[r0:r1].upper()
    p = _common_prefix_len(old, prefix.upper())
    q = _common_prefix_len(old[p:][::-1], suffix.upper()[::-1])
    a0, a1 = r0 + p, r1 - q
    mid = new_part[p:len(new_part) - q]
    spacer_start = r0 + len(prefix)

    record     = copy.deepcopy(frame)
    record.seq = Seq(seq[:a0] + mid + seq[a1:])
    record.features = [
        f for f in (_remap_feature(f, a0, a1, len(mid)) for f in frame.features)
        if f is not None
    ]
    record.features.append(
        SeqFeature(
            FeatureLocation(spacer_start, spacer_start + len(sgrna_seq), strand=+1),
            type="misc_feature",
            qualifiers={
                "label": [sgrna_id],
                "note":  [f"sgRNA sequence: {sgrna_seq}"],
                "gene":  [sgrna_id],
            },
        )
    )
    if reverse:
        record = record.reverse_complement(**_RC_KEEP)

    record.id          = sgrna_id
    record.name        = sgrna_id[:16]   # LOCUS name field: max 16 chars
    record.description = f"{sgrna_id} cloned into {vector_record.name.strip()}"
    record.features.sort(key=lambda f: int(f.location.start))
    return record


def build_plasmid_records(
    rows: list[dict],
    vector_record,
    site: CloningSite,
    seq_col: str = "sgRNA_seq",
    id_col: str  = "sgRNA_id",
) -> list:
    """
    Generate one annotated SeqRecord per row. Rows must already have been
    through :func:`add_cloning_oligos` (they carry ``insert_core``).
    """
    return [
        insert_sgrna(vector_record, site, row["insert_core"], row[seq_col], row[id_col])
        for row in rows
    ]


def named_record(record, vector_name: str):
    """
    Copy of a plasmid record named '{vector_name}_{sgRNA_id}' (LOCUS name and accession), so that
    SnapGene and Benchling show the full construct name, e.g. 'ipUSEPR_sg_Luc_474'.
    """
    r             = copy.deepcopy(record)
    full_id       = f"{sanitize_name(vector_name, 'vector')}_{record.id}"
    r.id          = full_id
    r.name        = full_id
    r.description = f"{record.id} cloned into {sanitize_name(vector_name, 'vector')}"
    return r


def genbank_string(records: list) -> str:
    """Serialise records as GenBank text (one or several records, in memory)."""
    import io
    from Bio import SeqIO

    buf = io.StringIO()
    SeqIO.write(records, buf, "genbank")
    return buf.getvalue()


def write_gb_files(
    records: list,
    output_dir: Path,
    vector_name: str,
    multi_gb: bool = False,
    combined_stem: str = "combined_plasmids",
) -> None:
    """
    Write GenBank files, one per sgRNA (default) or a single combined file.

    Each record's LOCUS name and accession are set to '{vector_name}_{sgRNA_id}'
    (see :func:`named_record`).
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    combined_stem = sanitize_name(combined_stem, "combined_plasmids")
    named = [named_record(r, vector_name) for r in records]

    if multi_gb:
        (output_dir / f"{combined_stem}.gb").write_text(genbank_string(named))
    else:
        for record in named:
            (output_dir / f"{record.id}.gb").write_text(genbank_string([record]))


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def read_input(path: Path) -> list[dict]:
    """
    Read an input spreadsheet or delimited file into a list of row dicts.
    Supports .xlsx / .xlsm (via openpyxl), .csv, and .tsv.

    Completely empty rows (all values empty or whitespace-only) and columns with a blank
    header are dropped; both are common in files saved from Excel. Raises ValueError for
    an unsupported type, an empty table, or more than ``MAX_SGRNAS`` rows.
    """
    suffix = path.suffix.lower()

    if suffix in {".xlsx", ".xlsm"}:
        import openpyxl
        wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
        try:
            rows_iter = wb.active.iter_rows(values_only=True)
            raw_headers = next(rows_iter, None)
            if raw_headers is None:
                raise ValueError("The spreadsheet is empty.")
            headers = [str(h).strip() if h is not None else "" for h in raw_headers]
            rows = [
                {h: (str(v) if v is not None else "") for h, v in zip(headers, row) if h}
                for row in rows_iter
            ]
        finally:
            wb.close()
    elif suffix in {".csv", ".tsv"}:
        with open(path, newline="", encoding="utf-8-sig") as fh:
            rows = [{k: v for k, v in row.items() if k}
                    for row in csv.DictReader(fh, delimiter="," if suffix == ".csv" else "\t")]
    else:
        raise ValueError(f"Unsupported file type: '{suffix}'. Use .xlsx, .csv, or .tsv")

    rows = [r for r in rows if any(str(v or "").strip() for v in r.values())]
    if not rows:
        raise ValueError("No sgRNA rows were found in the table.")
    if len(rows) > MAX_SGRNAS:
        raise ValueError(
            f"The table has {len(rows)} rows; this app handles up to {MAX_SGRNAS} sgRNAs at a time. "
            "Split the table and run it in batches."
        )
    return rows


def write_idt_excel(idt_rows: list[dict], path: Path) -> None:
    """
    Write IDT-formatted rows to Excel, replicating the official IDT bulk
    upload template (bold Calibri headers, IDT column widths, no index).
    """
    import openpyxl
    from openpyxl.styles import Font

    path.parent.mkdir(parents=True, exist_ok=True)

    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Sheet1"

    header_font = Font(name="Calibri", bold=True)
    for col_idx, col_name in enumerate(IDT_COLUMNS, start=1):
        ws.cell(row=1, column=col_idx, value=col_name).font = header_font

    data_font = Font(name="Calibri")
    for row_idx, row in enumerate(idt_rows, start=2):
        for col_idx, col_name in enumerate(IDT_COLUMNS, start=1):
            cell = ws.cell(row=row_idx, column=col_idx, value=row[col_name])
            cell.font = data_font
            if isinstance(row[col_name], str):
                cell.data_type = "s"      # a leading '=' must stay text, never a formula

    for col_letter, width in IDT_COL_WIDTHS.items():
        ws.column_dimensions[col_letter].width = width

    wb.save(path)

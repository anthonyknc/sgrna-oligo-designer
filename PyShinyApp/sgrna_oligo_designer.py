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
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# IDT bulk upload template (matches the IDT spec exactly)
IDT_COLUMNS    = ["Name", "Sequence", "Scale", "Purification"]
IDT_SCALE      = "25nm"
IDT_PURIF      = "STD"
IDT_COL_WIDTHS = {"A": 12.54, "B": 41.54, "C": 10.18, "D": 18.82}

# The enzymes listed first in the app (BsmBI/Esp3I are isoschizomers). Any other eligible
# enzyme can be chosen under "Other"; see all_enzymes(). Cut geometry always comes from
# Bio.Restriction.
ENZYMES = ("BsmBI", "Esp3I", "BsaI", "BbsI")

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

# Largest sgRNA table accepted in one run (keeps the web app responsive).
MAX_SGRNAS = 1000

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


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

    @property
    def removed_length(self) -> int:
        """Bases between the two overhangs that are replaced by the insert."""
        return self.replace_end - self.replace_start


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
        # the shorter piece is the one crossing the origin
        raise ValueError(
            f"The two {enzyme.name} sites split the plasmid into a {removed} bp piece between the "
            f"sites and a {other_piece} bp piece that crosses the origin (position 1) of the map. "
            "The shorter piece is taken to be the stuffer, and a stuffer that crosses the origin is "
            "not supported. Move the origin elsewhere in your plasmid editor, re-export the "
            "GenBank file, and try again. If that does not fit your vector, check that "
            f"{enzyme.name} is the enzyme used to clone into it (it may just happen to cut twice "
            "elsewhere), or that the stuffer is not the longer piece, which is not supported yet."
        )
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
    )


def _looks_like_scaffold_start(right_overhang: str) -> bool:
    """
    True if a right-hand overhang is the start of the Cas9 scaffold (GTTT...). Needs at
    least 3 bases to be informative; longer overhangs are compared on their first 4 bases.
    """
    n = min(len(right_overhang), len(SCAFFOLD_START))
    return len(right_overhang) >= 3 and right_overhang[:n] == SCAFFOLD_START[:n]


def find_cloning_site(vector_record, enzyme_name: str, orientation: str = "auto") -> CloningSite:
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

    Raises
    ------
    ValueError
        If the vector is not circular, does not have exactly two sites for the
        enzyme, or the cloning site spans the plasmid origin.
    """
    if orientation not in ORIENTATIONS:
        raise ValueError(f"orientation must be one of {ORIENTATIONS}, got '{orientation}'.")
    enzyme = get_enzyme(enzyme_name)
    if vector_record.annotations.get("topology") != "circular":
        raise ValueError(
            "The vector GenBank file is not marked circular. Only circular plasmids are supported."
        )

    seq    = str(vector_record.seq).upper()
    rc_seq = reverse_complement(seq)
    fx0, fy0 = _locate(seq, enzyme)                  # raises with a helpful message

    def forward(warnings=()):
        return _site_from(seq, enzyme, fx0, fy0, "forward", warnings)

    def reverse():
        rx0, ry0 = _locate(rc_seq, enzyme)
        return _site_from(rc_seq, enzyme, rx0, ry0, "reverse")

    if orientation == "forward":
        return forward()
    if orientation == "reverse":
        return reverse()

    fwd = forward()
    if _looks_like_scaffold_start(fwd.right_overhang):
        return fwd
    rev = reverse()
    if _looks_like_scaffold_start(rev.right_overhang):
        return rev
    return forward((
        f"Could not confirm the sgRNA orientation: expected the scaffold start ({SCAFFOLD_START}) "
        f"at the start of the right-hand overhang, found {fwd.right_overhang} (forward) / {rev.right_overhang} "
        "(reverse). Assuming forward; set the orientation under Advanced options if this is wrong.",
    ))


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
    Build the sense and antisense oligos for one insert.

    sense     = left overhang + core
    antisense = reverse complement of (core + right overhang)
    """
    return (
        site.left_overhang + core,
        reverse_complement(core + site.right_overhang),
    )


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
    site_pats = (site.enzyme_site, reverse_complement(site.enzyme_site))

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
        window = site.left_context + site.left_overhang + core + site.right_overhang + site.right_context
        if any(p in window for p in site_pats):
            notes.append(
                f"contains a {site.enzyme} recognition site ({site.enzyme_site}); the enzyme "
                "would cut this insert during cloning, so this guide will likely fail."
            )

        row[seq_col]           = seq
        row["insert_core"]     = core
        row["oligo_sense"], row["oligo_antisense"] = design_oligos(core, site)
        row["warnings"]        = [f"sgRNA '{sgrna_id}' {n}" for n in notes]
    return rows


_NAME_OK  = re.compile(r"[A-Za-z0-9_.\-]+")
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
            "letters, digits, underscores (_), hyphens (-) and dots (.), with no spaces. "
            "Rename them in your table, for example by replacing spaces with underscores."
        )
    duplicates = [i for i, n in Counter(r[id_col] for r in rows).items() if n > 1]
    if duplicates:
        raise ValueError(
            f"Duplicate sgRNA IDs found: {duplicates}. "
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
    # If the vector already has the 5' G right after the left overhang (as ipUSEPR and
    # lentiCRISPR do), keep that base and its annotations; replace only what follows.
    if g_len and seq[r0:r0 + g_len].upper() == core[:g_len]:
        r0 += g_len
        new_part = core[g_len:]
    else:
        new_part = core
    spacer_start = r0 if new_part == sgrna_seq else r0 + g_len

    record     = copy.deepcopy(frame)
    record.seq = Seq(seq[:r0] + new_part + seq[r1:])
    record.features = [
        f for f in (_remap_feature(f, r0, r1, len(new_part)) for f in frame.features)
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
    so that SnapGene and Benchling display the full construct name
    (e.g. 'ipUSEPR_sg_Luc_474').
    """
    from Bio import SeqIO

    output_dir.mkdir(parents=True, exist_ok=True)
    vector_name   = sanitize_name(vector_name, "vector")
    combined_stem = sanitize_name(combined_stem, "combined_plasmids")

    def _named(record):
        r             = copy.deepcopy(record)
        full_id       = f"{vector_name}_{record.id}"
        r.id          = full_id
        r.name        = full_id
        r.description = f"{record.id} cloned into {vector_name}"
        return r

    if multi_gb:
        out_path = output_dir / f"{combined_stem}.gb"
        with open(out_path, "w") as fh:
            SeqIO.write([_named(r) for r in records], fh, "genbank")
    else:
        for record in records:
            named = _named(record)
            with open(output_dir / f"{named.id}.gb", "w") as fh:
                SeqIO.write([named], fh, "genbank")


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

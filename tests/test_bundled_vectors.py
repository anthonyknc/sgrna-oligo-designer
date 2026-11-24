"""
The real Addgene vectors bundled with the app: each must be recognised (promoter, enzymes) and
give the oligo design published for that vector. The expected overhangs and flanks are the ones
in the vector's own protocol (Ran 2013 for px330, Sanjana 2014 for the lentiGuide/lentiCRISPR
family, the Weissman/Bassik BstXI + BlpI oligo design for the two-enzyme vectors).
"""
from pathlib import Path

import pytest

import sgrna_oligo_designer as core

VECTORS = Path(__file__).resolve().parent.parent / "PyShinyApp" / "vectors"
S = "TACGATTCGATACGCTAGCA"                     # a spacer without sites for any of these enzymes
RC = core.reverse_complement(S)

# file, enzymes, promoter, g_mode, expected sense oligo, expected antisense oligo
CASES = [
    ("ipUSEPR.gb",         "BsmBI", "human U6", "always", "CACCG" + S, "AAAC" + RC + "C"),
    ("lentiCRISPRv2.gb",   "BsmBI", "human U6", "always", "CACCG" + S, "AAAC" + RC + "C"),
    ("lentiGuide-Puro.gb", "BsmBI", "human U6", "always", "CACCG" + S, "AAAC" + RC + "C"),
    ("pX330.gb",           "BbsI",  "human U6", "always", "CACCG" + S, "AAAC" + RC + "C"),
    ("pU6-BbsI_CBh-Cas9-T2A-mCherry-H1.gb", "BbsI", "human U6", "always", "CACCG" + S, "AAAC" + RC + "C"),
    ("phH1-gRNA.gb",       "BbsI",  "H1",       "never",  "TCCC" + S, "AAAC" + RC),
    ("pU6-sgRNA_EF1Alpha-puro-T2A-BFP.gb", "BstXI + BlpI", "mouse U6", "never",
     "TTG" + S + "GTTTAAGAGC", "TTAGCTCTTAAAC" + RC + "CAACAAG"),
    ("pMCB306.gb",         "BstXI + BlpI", "mouse U6", "never",
     "TTG" + S + "GTTTAAGAGC", "TTAGCTCTTAAAC" + RC + "CAACAAG"),
]


@pytest.mark.parametrize("file, enzymes, promoter, g_mode, sense, antisense", CASES)
def test_vector(file, enzymes, promoter, g_mode, sense, antisense):
    rec = core.read_vector_gb(VECTORS / file)
    best = core.suggest_cloning_sites(rec)[0]
    assert best.enzymes == enzymes and best.promoter.startswith(promoter)

    one = best.enzyme_b == ""
    site = core.find_cloning_site(rec, best.enzyme, **({} if one else {"enzyme_b": best.enzyme_b}))
    rows = [{"sgRNA_id": "t1", "sgRNA_seq": S}]
    core.add_cloning_oligos(rows, site, g_mode=g_mode)
    assert (rows[0]["oligo_sense"], rows[0]["oligo_antisense"]) == (sense, antisense)
    assert rows[0]["warnings"] == []

    (plasmid,) = core.build_plasmid_records(rows, rec, site)
    seq = str(plasmid.seq)
    assert seq.count(S) == 1
    assert seq[seq.index(S) + len(S):].startswith(("GTTT", "GTTTAAGAGC"))     # the scaffold follows


def test_every_bundled_vector_is_listed_in_the_app():
    import app
    listed = {path.name for _label, path in app._VECTORS.values()}
    assert listed == {p.name for p in VECTORS.glob("*.gb")}

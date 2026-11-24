"""Pol III promoter detection and promoter-anchored suggestions."""
import random
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import sgrna_oligo_designer as core
from promoters import find_pol3_promoters, _from_sequence

ROOT = Path(__file__).resolve().parent.parent
VECTORS = ROOT / "PyShinyApp" / "vectors"
PMCB306 = ROOT / "tests" / "fixtures" / "pMCB306.gb"


def load(path):
    return core.read_vector_gb(Path(path))


def bare(record):
    """The same plasmid without any annotations."""
    r = SeqRecord(record.seq, id=record.id, name=record.name)
    r.annotations["topology"] = "circular"
    r.annotations["molecule_type"] = "DNA"
    return r


@pytest.mark.parametrize("name, strand", [
    ("ipUSEPR.gb", +1), ("lentiCRISPRv2.gb", +1), ("LentiCRISPRv2GFP.gb", -1)])
def test_annotated_human_u6(name, strand):
    (p,) = find_pol3_promoters(load(VECTORS / name))
    assert (p.name, p.strand, p.source) == ("human U6", strand, "annotation")


@pytest.mark.parametrize("name", ["ipUSEPR.gb", "lentiCRISPRv2.gb", "LentiCRISPRv2GFP.gb"])
def test_sequence_fallback_matches_annotation(name):
    rec = load(VECTORS / name)
    (annotated,) = find_pol3_promoters(rec)
    (found,) = find_pol3_promoters(bare(rec))
    assert found.source == "sequence" and found.name == "human U6" and found.strand == annotated.strand
    # same 3' end as the annotation, give or take the CACC the reference leaves out
    assert abs(found.in_frame(len(rec.seq), "forward" if found.strand == 1 else "reverse")
               - annotated.in_frame(len(rec.seq), "forward" if annotated.strand == 1 else "reverse")) <= 5


def test_mouse_u6_found_by_sequence_only_in_pmcb306():
    (p,) = find_pol3_promoters(load(PMCB306))
    assert (p.name, p.source, p.strand) == ("mouse U6", "sequence", +1)


def test_no_cross_species_or_random_hits():
    random.seed(1)
    junk = "".join(random.choice("ACGT") for _ in range(8000))
    assert _from_sequence(junk) == []
    assert [p.name for p in _from_sequence(str(load(VECTORS / "ipUSEPR.gb").seq))] == ["human U6"]


@pytest.mark.parametrize("label, name", [
    ("hU6 promoter", "human U6"), ("mouse U6 promoter", "mouse U6"), ("U6", "U6"),
    ("H1 promoter", "H1"), ("7SK promoter", "7SK"), ("Pol III promoter", "Pol III promoter")])
def test_annotation_names(label, name):
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    rec = bare(load(VECTORS / "ipUSEPR.gb"))
    rec.features.append(SeqFeature(FeatureLocation(10, 200, strand=+1), type="promoter",
                                   qualifiers={"label": [label]}))
    assert [p.name for p in find_pol3_promoters(rec)] == [name]


def test_primers_and_other_features_are_not_promoters():
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    rec = bare(load(VECTORS / "ipUSEPR.gb"))
    rec.features.append(SeqFeature(FeatureLocation(10, 30, strand=+1), type="primer_bind",
                                   qualifiers={"label": ["hU6-F"]}))
    rec.features.append(SeqFeature(FeatureLocation(40, 90, strand=+1), type="CDS",
                                   qualifiers={"label": ["U6 binding protein"]}))
    # no annotation matched, so detection falls back to the sequence (which does find the real U6)
    assert [p.source for p in find_pol3_promoters(rec)] == ["sequence"]


@pytest.mark.parametrize("path, enzyme, kw", [
    (VECTORS / "ipUSEPR.gb", "BsmBI", {}),
    (VECTORS / "LentiCRISPRv2GFP.gb", "BsmBI", {}),
    (PMCB306, "BstXI", {"enzyme_b": "BlpI"})])
def test_site_records_its_promoter(path, enzyme, kw):
    site = core.find_cloning_site(load(path), enzyme, **kw)
    assert "U6 promoter" in site.promoter


def test_suggestions():
    first = lambda path: [s.enzymes for s in core.suggest_cloning_sites(load(path))]
    assert first(VECTORS / "ipUSEPR.gb")[0] == "BsmBI"
    assert first(VECTORS / "LentiCRISPRv2GFP.gb")[0] == "BsmBI"
    assert first(PMCB306) == ["BstXI + BlpI"]


def test_no_promoter_no_suggestions():
    random.seed(2)
    rec = bare(load(VECTORS / "ipUSEPR.gb"))
    rec.seq = Seq("".join(random.choice("ACGT") for _ in range(3000)))
    assert core.suggest_cloning_sites(rec) == []


def test_promoter_decides_orientation_when_scaffold_cannot():
    rec = load(VECTORS / "LentiCRISPRv2GFP.gb")          # sgRNA on the reverse strand
    site = core.find_cloning_site(rec, "BsmBI", "reverse")
    assert site.right_overhang == "GTTT"
    seq, n = str(rec.seq), len(rec.seq)
    j = site.replace_end                                  # right overhang, sense frame
    mutated = seq[:n - j - 4] + core.reverse_complement("GCTT") + seq[n - j:]
    assert core.find_cloning_site(SeqRecord(Seq(mutated), annotations=dict(rec.annotations)),
                                  "BsmBI", "reverse").right_overhang == "GCTT"
    rec2 = SeqRecord(Seq(mutated), id=rec.id, name=rec.name, features=rec.features,
                     annotations=dict(rec.annotations))
    auto = core.find_cloning_site(rec2, "BsmBI", "auto")
    assert auto.orientation == "reverse"
    assert not any("Could not confirm" in w for w in auto.warnings)

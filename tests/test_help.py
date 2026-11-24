"""The help/*.md files are the single source of truth for the in-app Guide and the README."""

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
HELP = ROOT / "PyShinyApp" / "help"
NAMES = ["getting_started", "how_it_works", "troubleshooting", "citing"]


@pytest.mark.parametrize("name", NAMES)
def test_help_file_has_content(name):
    text = (HELP / f"{name}.md").read_text(encoding="utf-8")
    assert text.startswith("# ") and len(text) > 200


@pytest.mark.parametrize("name", NAMES)
def test_help_links_are_absolute(name):
    """Relative links work on GitHub but break inside the app."""
    text = (HELP / f"{name}.md").read_text(encoding="utf-8")
    bad = re.findall(r"\]\((?!https?://)[^)]*\)", text)
    assert not bad, bad


def test_help_files_match_the_current_design():
    """Guard against stale text from the old stuffer/CLI design."""
    text = "\n".join((HELP / f"{n}.md").read_text(encoding="utf-8") for n in NAMES)
    for stale in ("--stuffer", "--vector", "pixi", "python3 sgrna", "sense prefix", "Stuffer sequence /"):
        assert stale not in text, stale


def test_readme_links_every_help_file_and_repeats_the_tool_citation():
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    for n in NAMES:
        assert f"PyShinyApp/help/{n}.md" in readme
    citation = next(l for l in (HELP / "citing.md").read_text(encoding="utf-8").splitlines()
                    if l.startswith("> Chan, A.K.N."))
    assert citation in readme


def test_app_renders_every_guide_page():
    import app
    for n in NAMES:
        assert "<h" in str(app._guide(n))


def test_app_version_matches_pyproject_and_last_update_is_a_date():
    import datetime
    import tomllib

    import app
    version = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]["version"]
    assert app._VERSION == f"v{version}"
    datetime.date.fromisoformat(app._LAST_UPDATE)      # raises unless YYYY-MM-DD


def test_guide_states_the_main_points_correctly():
    how = (HELP / "how_it_works.md").read_text(encoding="utf-8")
    assert "strongly prefers" in how and "requires a" not in how          # 5' G wording
    assert "Do the oligos need to be phosphorylated?" in how               # PNK note
    assert "No. The oligos are ordered without a 5' phosphate" in how

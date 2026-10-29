"""Release wiring checks; scientific verification lives in numerical tests."""

from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]


def test_pages_workflow_tests_before_deploying_committed_site():
    path = ROOT / ".github/workflows/pages.yml"
    assert path.is_file(), "Pages workflow has not been implemented"
    workflow = path.read_text()
    for command in (
        "uv sync --frozen",
        "uv run pytest -q",
        "node --test tests/test_site.mjs",
        "uv run ruff check src tests scripts",
        "uv run ruff format --check src tests scripts",
    ):
        assert command in workflow
    assert "needs: verify" in workflow
    assert "path: site" in workflow
    actions = re.findall(r"uses:\s+(\S+)", workflow)
    assert actions
    assert all(re.fullmatch(r"[\w/-]+@[0-9a-f]{40}", action) for action in actions)
    assert "contents: read" in workflow
    assert "pages: write" in workflow
    assert "id-token: write" in workflow
    assert "pull_request_target" not in workflow


def test_public_docs_and_shell_have_relative_reproduction_paths():
    for name in ("README.md", "docs/PASSIVE.md", "docs/ACTIVE.md", "docs/PIPELINE.md"):
        text = (ROOT / name).read_text()
        assert "/home/nish/" not in text, name
    assert 'href="./styles.css"' in (ROOT / "site/index.html").read_text()
    assert "data/reference/" in (ROOT / ".gitignore").read_text()


def test_default_uv_python_matches_the_accepted_ci_interpreter():
    version = ROOT / ".python-version"
    assert version.is_file(), "Default reproduction interpreter is not pinned"
    assert version.read_text().strip() == "3.11"

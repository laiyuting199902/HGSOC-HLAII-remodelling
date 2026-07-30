from pathlib import Path
import csv
import re


ROOT = Path(__file__).resolve().parents[1]


def iter_text_files():
    suffixes = {".py", ".R", ".md", ".tsv", ".csv", ".json", ".txt", ".yml", ".yaml", ".sh", ".cpp"}
    for path in ROOT.rglob("*"):
        if ".git" in path.parts or not path.is_file() or path.suffix not in suffixes:
            continue
        yield path


def test_repository_contains_code_and_supporting_metadata_only():
    forbidden_directories = {
        "source_data",
        "processed_summaries",
        "figures",
        "manuscript",
        "manuscripts",
        "submission",
        "submission_files",
        "publication_figures",
    }
    assert not forbidden_directories.intersection(path.name for path in ROOT.iterdir() if path.is_dir())

    forbidden_suffixes = {".doc", ".docx", ".pdf", ".png", ".jpg", ".jpeg", ".svg", ".tif", ".tiff"}
    artifacts = [
        path.relative_to(ROOT)
        for path in ROOT.rglob("*")
        if ".git" not in path.parts and path.is_file() and path.suffix.lower() in forbidden_suffixes
    ]
    assert not artifacts, f"generated or manuscript artifacts remain: {artifacts}"


def test_no_local_absolute_paths_or_repository_placeholders():
    forbidden = ("/" + "Users/", "/" + "Volumes/")
    for path in iter_text_files():
        content = path.read_text(encoding="utf-8", errors="ignore")
        for needle in forbidden:
            assert needle not in content, f"local path remains in {path}: {needle}"
        assert ("ZEN" + "ODO DOI") not in content.upper()


def test_no_obvious_embedded_credentials():
    pattern = re.compile(r"(?i)(api[_-]?key|access[_-]?token|password)\s*[:=]\s*['\"][^'\"]+['\"]")
    for path in iter_text_files():
        assert pattern.search(path.read_text(encoding="utf-8", errors="ignore")) is None


def test_no_document_editing_dependency_or_script():
    requirements = (ROOT / "environment" / "requirements.txt").read_text(encoding="utf-8").lower()
    assert "python-docx" not in requirements
    for path in (ROOT / "scripts").iterdir():
        assert "submission" not in path.name.lower()
        assert "manuscript" not in path.name.lower()
        if path.suffix == ".py":
            content = path.read_text(encoding="utf-8", errors="ignore")
            assert "from docx" not in content
            assert "import docx" not in content


def test_workflow_manifest_references_existing_scripts():
    with (ROOT / "manifests" / "analysis_workflow.tsv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows
    for row in rows:
        assert (ROOT / "scripts" / row["script"]).is_file(), row["script"]

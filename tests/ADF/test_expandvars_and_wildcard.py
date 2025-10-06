import pathlib as pl

from pyfrag.executables.adf.input import expandvars_backslash
from pyfrag.executables.adf.mol_handling import find_files_with_wildcard_option


def test_expandvars_backslash_env(monkeypatch):
    # Mock environment variable
    monkeypatch.setenv("SLURM_SUBMIT_DIR", "/mock/dir")
    # Path with env variable
    path = pl.Path("$SLURM_SUBMIT_DIR/testfile.xyz")
    # Should expand to /mock/dir/testfile.xyz
    expanded = expandvars_backslash(path)
    assert str(expanded) == "/mock/dir/testfile.xyz"


def test_expandvars_backslash_escaped(monkeypatch):
    monkeypatch.setenv("SLURM_SUBMIT_DIR", "/mock/dir")
    # Escaped $ should not expand
    path = r"\$SLURM_SUBMIT_DIR/testfile.xyz"
    expanded = expandvars_backslash(path)
    # The backslash is preserved, and $SLURM_SUBMIT_DIR is not expanded
    assert str(expanded) == r"\$SLURM_SUBMIT_DIR/testfile.xyz"


def test_find_files_with_wildcard_option(tmp_path):
    # Create files
    file1 = tmp_path / "coord1.xyz"
    file2 = tmp_path / "coord2.xyz"
    file1.write_text("test1")
    file2.write_text("test2")
    # Use wildcard
    pattern = tmp_path / "coord*.xyz"
    result = find_files_with_wildcard_option([pattern])
    result_names = sorted([f.name for f in result])
    assert result_names == ["coord1.xyz", "coord2.xyz"]


def test_find_files_with_wildcard_option_no_wildcard(tmp_path):
    file1 = tmp_path / "coord1.xyz"
    file1.write_text("test1")
    result = find_files_with_wildcard_option([file1])
    assert result == [file1]

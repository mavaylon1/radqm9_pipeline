from pathlib import Path

from ase.io import read

import pytest


@pytest.fixture(scope="session")
def test_dir():
    return Path(__file__).parent.parent.joinpath("test_files").resolve()

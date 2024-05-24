import json
import datetime

import pytest

import numpy as np

from monty.io import zopen

from emmet.core.qchem.task import TaskDocument

from radqm9_pipeline.models.force import ForcePointDoc


try:
    from openbabel.openbabel import OBAlign

    _ = OBAlign()
    has_eigen = True
except ImportError:
    has_eigen = False


@pytest.fixture(scope="session")
def force_calc(test_dir):
    with zopen(test_dir / "force_calc.json.gz") as f:
        data = json.load(f)
    
    data["last_updated"] = datetime.datetime.now()

    task = TaskDocument(**data)
    return task


def test_force_point(force_calc):
    doc = ForcePointDoc.from_task(
        force_calc, molecule_id="b9ba54febc77d2a9177accf4605767db-C1Li2O3-1-2"
    )
    
    # Basic properties
    assert doc.property_name == "force_single_point"
    assert doc.species == ['O', 'C', 'O', 'C', 'C', 'O', 'C', 'O', 'O', 'H', 'H']
    
    # Forces
    assert np.allclose(np.array(doc.forces), np.array(doc.precise_forces), atol=0.001)
    assert doc.cds_forces is None

    # Partial charges/spins
    assert np.allclose(
        np.array(doc.mulliken_partial_charges),
        np.array(
            [
                -0.536376,
                0.454763,
                0.235882,
                -1.183636,
                -1.18708,
                0.443149,
                -0.429893,
                0.440251,
                0.237663,
                0.276062,
                0.249216
            ]
        )
    )
    assert np.allclose(
        np.array(doc.nbo_partial_spins),
        np.array(
            [
                0.0014,
                0.00849,
                0.05366,
                0.3188,
                0.31882,
                0.10078,
                0.03775,
                0.10064,
                0.05373,
                0.00225,
                0.00367
            ]
        )
    )

    # Dipole moment
    assert doc.dipole_moment == (pytest.approx(0.1388), pytest.approx(-12.5833), pytest.approx(-1.2874))


def test_schema():
    ForcePointDoc.model_json_schema()

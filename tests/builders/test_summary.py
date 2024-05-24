import pytest

from maggma.stores import JSONStore, MemoryStore

from emmet.builders.qchem.molecules import MoleculesAssociationBuilder, MoleculesBuilder
from emmet.builders.molecules.atomic import PartialChargesBuilder, PartialSpinsBuilder
from emmet.builders.molecules.electric import ElectricMultipoleBuilder
from emmet.builders.molecules.redox import RedoxBuilder
from emmet.builders.molecules.thermo import ThermoBuilder
from emmet.builders.molecules.vibration import VibrationBuilder

from radqm9_pipeline.builders.summary import SummaryBuilder


@pytest.fixture(scope="session")
def tasks(test_dir):
    return JSONStore(test_dir / "force_traj_tasks.json.gz")


@pytest.fixture(scope="session")
def mols(tasks):
    assoc_store = MemoryStore(key="molecule_id")
    stage_one = MoleculesAssociationBuilder(tasks=tasks, assoc=assoc_store)
    stage_one.run()

    return assoc_store


@pytest.fixture(scope="session")
def charges(test_dir):
    return MemoryStore(key="molecule_id")


@pytest.fixture(scope="session")
def spins(test_dir):
    return MemoryStore(key="molecule_id")


@pytest.fixture(scope="session")
def multipoles(test_dir):
    return MemoryStore(key="molecule_id")


@pytest.fixture(scope="session")
def redox(test_dir):
    return MemoryStore(key="molecule_id")


@pytest.fixture(scope="session")
def thermo(test_dir):
    return MemoryStore(key="molecule_id")


@pytest.fixture(scope="session")
def vibes(test_dir):
    return MemoryStore(key="molecule_id")


@pytest.fixture(scope="session")
def summary():
    return MemoryStore(key="molecule_id")


def test_summary_one(
    tasks,
    mols,
    charges,
    spins,
    multipoles,
    redox,
    thermo,
    vibes,
    summary,
):

    charge_build = PartialChargesBuilder(tasks, mols, charges)
    charge_build.run()

    spins_build = PartialSpinsBuilder(tasks, mols, spins)
    spins_build.run()

    multipole_build = ElectricMultipoleBuilder(tasks, mols, multipoles)
    multipole_build.run()

    thermo_build = ThermoBuilder(tasks, mols, thermo)
    thermo_build.run()

    redox_build = RedoxBuilder(tasks, mols, thermo, redox)
    redox_build.run()

    vibe_build = VibrationBuilder(tasks, mols, vibes)
    vibe_build.run()

    builder = SummaryBuilder(
        molecules=mols,
        charges=charges,
        spins=spins,
        multipoles=multipoles,
        redox=redox,
        thermo=thermo,
        vibes=vibes,
        summary=summary,
    )
    builder.run()

    assert summary.count() == 8

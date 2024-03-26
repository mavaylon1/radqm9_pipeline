import pytest
from maggma.stores import JSONStore, MemoryStore

from emmet.builders.qchem.molecules import MoleculesAssociationBuilder, MoleculesBuilder

from radqm9_pipeline.builders.trajectory import TrajectoryBuilder


__author__ = "Evan Spotte-Smith <ewcspottesmith@lbl.gov>"


@pytest.fixture(scope="session")
def tasks_store(test_dir):
    return JSONStore(test_dir / "force_traj_tasks.json.gz")


@pytest.fixture(scope="session")
def mol_store(tasks_store):
    assoc_store = MemoryStore(key="molecule_id")
    stage_one = MoleculesAssociationBuilder(tasks=tasks_store, assoc=assoc_store)
    stage_one.run()

    return assoc_store


@pytest.fixture(scope="session")
def traj_store():
    return MemoryStore()


def test_trajectory_builder(tasks_store, mol_store, traj_store):
    builder = TrajectoryBuilder(tasks_store, mol_store, traj_store)
    builder.run()

    assert traj_store.count() == 8

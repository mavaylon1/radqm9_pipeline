import pytest
from maggma.stores import JSONStore, MemoryStore

from emmet.builders.qchem.molecules import MoleculesAssociationBuilder, MoleculesBuilder

from radqm9_pipeline.build.force import ForcesBuilder


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
def force_store():
    return MemoryStore()


def test_force_builder(tasks_store, mol_store, force_store):
    builder = ForcesBuilder(tasks_store, mol_store, force_store)
    builder.run()

    assert force_store.count() == 48
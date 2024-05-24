import pytest
from monty.serialization import loadfn

from radqm9_pipeline.models.summary import RadQM9SummaryDoc


@pytest.fixture(scope="session")
def docs_data(test_dir):
    raw = loadfn(test_dir / "summary_input.json.gz")
    return raw


def test_summary_doc(docs_data):
    summary_doc = RadQM9SummaryDoc.from_docs(
        molecule_id="17bde69a88c35896bd48bf001a33cbc6-C6H7N1O1-0-1", docs=docs_data
    )

    assert summary_doc.property_name == "summary"
    assert summary_doc.electronic_energy is not None
    assert summary_doc.total_enthalpy is not None
    assert summary_doc.frequencies is not None
    assert summary_doc.electron_affinity is None

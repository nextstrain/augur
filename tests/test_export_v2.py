import pytest
from Bio import Phylo

from augur import export_v2


ZIKA_TREE = "(1_0087_PF:0.00021787);"

@pytest.fixture
def zika_tree(tmp_path):
    with open(str(tmp_path / 'tree.nwk'), 'w+') as f:
        f.write(ZIKA_TREE)
    return Phylo.read(str(tmp_path / 'tree.nwk'), 'newick')

@pytest.fixture
def zika_metadata():
    return ({'1_0087_PF': {'strain': '1_0087_PF', 'virus': 'zika', 'accession': 'KX447509', 'date': '2013-12-XX'}}, [])


@pytest.fixture
def zika_node_data():
    return {'nodes': {'1_0087_PF': {'branch_length': 0.0002178684358946914, 'clock_length': 0.0002178684358946914, 'date': '2013-12-23' }}}

@pytest.fixture
def zika_node_data_no_raise():
    return {'nodes': {'1_0087_PF': {'branch_length': 0.0002178684358946914, 'clock_length': 0.0002178684358946914}}}


def test_parse_node_data_and_metadata(monkeypatch, zika_tree, zika_node_data, zika_metadata):
    monkeypatch.setattr(export_v2, 'read_node_data', lambda x: zika_node_data)
    monkeypatch.setattr(export_v2, 'read_metadata', lambda x: zika_metadata)
    with pytest.raises(KeyError) as override_exception:
        _, node_attrs, _, _ = export_v2.parse_node_data_and_metadata(zika_tree, None, None)
    assert '1_0087_PF' in str(override_exception.value)

def test_parse_node_data_and_metadata_no_raise(monkeypatch, zika_tree, zika_node_data_no_raise, zika_metadata):
    monkeypatch.setattr(export_v2, 'read_node_data', lambda x: zika_node_data_no_raise)
    monkeypatch.setattr(export_v2, 'read_metadata', lambda x: zika_metadata)
    _, node_attrs, _, _ = export_v2.parse_node_data_and_metadata(zika_tree, None, None)
    test_strain = '1_0087_PF'
    assert node_attrs[test_strain]['accession'] == zika_metadata[0][test_strain]['accession']

#!/usr/bin/env python3
import pytest

from augur.clades import read_in_clade_definitions

def test_read_in_clade_definitions_simple():
    clades = read_in_clade_definitions("tests/data/clades/simple_clades.tsv")
    assert clades == {
        'Clade_1': [('ctpE', 80, 'D')],
        'Clade_2': [('nuc', 30641, 'T')],
        'Clade_3': [('nuc', 444295, 'A'), ('pks8', 633, 'T')]
    }

def test_read_in_clade_definitions_with_empty_lines():
    clades = read_in_clade_definitions("tests/data/clades/empty_lines_clades.tsv")
    assert clades == {
        'Clade_1': [('ctpE', 80, 'D')],
        'Clade_2': [('nuc', 30641, 'T')],
        'Clade_3': [('nuc', 444295, 'A'), ('pks8', 633, 'T')]
    }

def test_read_in_clade_definitions_with_comments():
    clades = read_in_clade_definitions("tests/data/clades/commented_clades.tsv")
    assert clades == {
        'Clade_1': [('ctpE', 80, 'D')],
        'Clade_2': [('nuc', 30641, 'T')],
        'Clade_3': [('nuc', 444295, 'A'), ('pks8', 633, 'T')]
    }

def test_read_in_clade_definitions_inherit_simple():
    clades = read_in_clade_definitions("tests/data/clades/inherit_clades.tsv")
    assert clades == {
        'Clade_1': [('ctpE', 80, 'D')],
        'Clade_2': [('nuc', 30641, 'T')],
        'Clade_3': [('nuc', 30641, 'T'), ('pks8', 633, 'T')]
    }

def test_read_in_clade_definitions_inherit_chained():
    clades = read_in_clade_definitions("tests/data/clades/inherit_chained_clades.tsv")
    assert clades == {
        'Clade_1': [('ctpE', 80, 'D')],
        'Clade_2': [('ctpE', 80, 'D'),('nuc', 30641, 'T')],
        'Clade_3': [('ctpE', 80, 'D'),('nuc', 30641, 'T'), ('pks8', 633, 'T')]
    }

def test_read_in_clade_definitions_inherit_cycle_error():
    with pytest.raises(ValueError):
        read_in_clade_definitions("tests/data/clades/inherit_cycle_clades.tsv")

def test_read_in_clade_definitions_multiple_inheritance_error():
    with pytest.raises(ValueError):
        read_in_clade_definitions("tests/data/clades/multiple_inheritance_clades.tsv")

def test_read_in_clade_definitions_inheritance_from_nonexistent_clade_error():
    with pytest.raises(ValueError):
        read_in_clade_definitions("tests/data/clades/nonexistent_clade_inheritance_clades.tsv")

def test_read_in_clade_definitions_inheritance_from_self_error():
    with pytest.raises(ValueError):
        read_in_clade_definitions("tests/data/clades/self_inherit_clades.tsv")

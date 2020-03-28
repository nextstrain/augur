"""Tests for clades assignment to nodes in a tree"""
import pytest

from augur.clades import (
    assign_clades,
    get_reference_sequence_from_root_node,
    is_node_in_clade,
    read_in_clade_definitions,
    register_arguments,
    run
)


def test_assign_clades():
    pass


def test_get_reference_sequence():
    pass


def test_is_node_in_clade():
    pass


def test_read_in_clade_definitions_tsv_valid():
    # TODO need a real world valid file
    result = read_in_clade_definitions('tests/data/clades.tsv')
    assert type(result) is dict


def test_read_in_clade_definitions_tsv_no_gene():
    with pytest.raises(AttributeError):
        read_in_clade_definitions('tests/data/titer_model/h3n2_titers_subset.tsv')


def test_read_in_clade_definitions_not_tsv_no_gene():
    with pytest.raises(AttributeError):
        read_in_clade_definitions('tests/data/distance_map_weight_per_site.json')


def test_read_in_clade_definition_file_not_found():
    with pytest.raises(FileNotFoundError):
        read_in_clade_definitions('hello.txt')


def test_register_arguments_for_command_options():
    pass


def test_run_clades():
    pass

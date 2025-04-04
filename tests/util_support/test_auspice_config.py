from augur.util_support.auspice_config import merge_configs
from augur.errors import AugurError
from augur.types import ValidationMode

import pytest

SKIP = ValidationMode.SKIP

class TestScalarMergeTitle:

    def test_overwrite(self):
        configs = [
            {'title': 'FIRST',},
            {'title': 'SECOND',},
        ]
        merged = merge_configs(configs, SKIP)
        assert merged['title'] == 'SECOND'

    def test_a_or_b(self):
        merged = merge_configs([{'title': 'FIRST',}, {}], SKIP)
        assert merged['title'] == 'FIRST'
        merged = merge_configs([{}, {'title': 'SECOND',}], SKIP)
        assert merged['title'] == 'SECOND'



class TestListMergeColorings:

    def test_simple_extend(self):
        configs = [
            {
                'colorings': [
                    {'key': 'a', 'title': 'A', 'type': 'categorical'},
                ]
            },
            {
                'colorings': [
                    {'key': 'b', 'title': 'B', 'type': 'categorical'},
                    {'key': 'c', 'title': 'C', 'type': 'categorical'},
                ]
            }
        ]
        merged = merge_configs(configs, SKIP)
        assert len(merged['colorings']) == 3
        assert [c['key'] for c in merged['colorings']] == ['a', 'b', 'c']

    def test_a_or_b(self):
        config = {
            'colorings': [
                {'key': 'a', 'title': 'A', 'type': 'categorical'},
            ]
        }
        merged1 = merge_configs([config, {}], SKIP)
        assert len(merged1['colorings']) == 1
        assert [c['key'] for c in merged1['colorings']] == ['a']

        merged2 = merge_configs([{}, config], SKIP)
        assert merged1 == merged2

    def test_overwrite(self):
        configs = [
            {
                'colorings': [
                    {'key': 'a', 'title': 'A', 'type': 'categorical'},
                ]
            },
            {
                'colorings': [
                    {'key': 'b', 'title': 'B', 'type': 'categorical'},
                    {'key': 'a', 'title': 'A_new', 'type': 'temporal'},
                ]
            }
        ]
        merged = merge_configs(configs, SKIP)
        assert len(merged['colorings']) == 2
        assert [c['key'] for c in merged['colorings']] == ['a', 'b']
        assert merged['colorings'][0] == configs[1]['colorings'][1]


class TestListMergeMetadataColumns:

    def test_simple_extend(self):
        configs = [
            {
                'metadata_columns': ['one', 'two', 'three']
            },
            {
                'metadata_columns': ['four', 'five', 'six']
            }
        ]
        merged = merge_configs(configs, SKIP)
        assert len(merged['metadata_columns']) == 6
        assert merged['metadata_columns'][3] == 'four'

    def test_overwrite(self):
        configs = [
            {
                'metadata_columns': ['one', 'two', 'three']
            },
            {
                'metadata_columns': ['four', 'two', 'one']
            }
        ]
        merged = merge_configs(configs, SKIP)
        assert len(merged['metadata_columns']) == 4
        assert merged['metadata_columns'] == ['one', 'two', 'three', 'four']


class TestDictMergeDisplayDefaults:

    def test_simple_extend(self):
        configs = [
            {
                'display_defaults': {
                    "map_triplicate": True
                }
            },
            {
                'display_defaults': {
                    "color_by": "country"
                }
            }
        ]
        merged = merge_configs(configs, SKIP)
        assert list(merged['display_defaults'].items()) == [('map_triplicate', True), ('color_by', 'country')] # order stable


    def test_a_or_b(self):
        config = {
            'display_defaults': {
                "map_triplicate": False
            }
        }
        merged1 = merge_configs([config, {}], SKIP)
        assert list(merged1['display_defaults'].items()) == [('map_triplicate', False)]

        merged2 = merge_configs([{}, config], SKIP)
        assert merged1==merged2


    def test_overwrite(self):
        configs = [
            {
                'display_defaults': {
                    "panels": ['map', 'frequencies'],
                    "map_triplicate": True,
                }
            },
            {
                'display_defaults': {
                    "color_by": "country",
                    "map_triplicate": False,
                    "panels": ['tree'],
                }
            }
        ]
        merged = merge_configs(configs, SKIP)
        assert list(merged['display_defaults'].keys()) == ['panels', 'map_triplicate', 'color_by']
        assert merged['display_defaults']['map_triplicate'] is False
        assert merged['display_defaults']['panels'] == ['tree'] # original list clobbered

    def test_non_dicts_are_error(self):
        configs = [
            {
                'display_defaults': {
                    "map_triplicate": True
                }
            },
            {
                'display_defaults': 'country'
            }
        ]
        with pytest.raises(AugurError):
            merged = merge_configs(configs, SKIP)
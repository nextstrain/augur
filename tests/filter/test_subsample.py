import pytest
import pandas as pd

import augur.filter.subsample
from augur.errors import AugurError


@pytest.fixture
def valid_metadata() -> pd.DataFrame:
    columns = ['strain', 'date', 'country']
    data = [
        ("SEQ_1","2020-01-XX","A"),
        ("SEQ_2","2020-02-01","A"),
        ("SEQ_3","2020-03-01","B"),
        ("SEQ_4","2020-04-01","B"),
        ("SEQ_5","2020-05-01","B")
    ]
    return pd.DataFrame.from_records(data, columns=columns).set_index('strain')


class TestSequencesPerGroup:
    @pytest.mark.parametrize(
        "target_max_value, counts_per_group, expected_sequences_per_group",
        [
            (3, [2, 2], 1),
            (3, [2, 1], 3),
            (9, [5, 5], 3),
            (9, [5, 4], 9),
            (9, [5, 3], 9),
        ],
    )
    def test_sequences_per_group(self, target_max_value, counts_per_group, expected_sequences_per_group):
        assert augur.filter.subsample._calculate_sequences_per_group(target_max_value, counts_per_group) == expected_sequences_per_group


class TestFilterGroupBy:
    def test_filter_groupby_strain_subset(self, valid_metadata: pd.DataFrame):
        metadata = valid_metadata.copy()
        strains = ['SEQ_1', 'SEQ_3', 'SEQ_5']
        group_by_strain = augur.filter.subsample.get_groups_for_subsampling(strains, metadata)
        assert group_by_strain == {
            'SEQ_1': ('_dummy',),
            'SEQ_3': ('_dummy',),
            'SEQ_5': ('_dummy',)
        }

    def test_filter_groupby_dummy(self, valid_metadata: pd.DataFrame):
        metadata = valid_metadata.copy()
        strains = metadata.index.tolist()
        group_by_strain = augur.filter.subsample.get_groups_for_subsampling(strains, metadata)
        assert group_by_strain == {
            'SEQ_1': ('_dummy',),
            'SEQ_2': ('_dummy',),
            'SEQ_3': ('_dummy',),
            'SEQ_4': ('_dummy',),
            'SEQ_5': ('_dummy',)
        }

    def test_filter_groupby_invalid_error(self, valid_metadata: pd.DataFrame):
        groups = ['invalid']
        metadata = valid_metadata.copy()
        strains = metadata.index.tolist()
        with pytest.raises(AugurError) as e_info:
            augur.filter.subsample.get_groups_for_subsampling(strains, metadata, group_by=groups)
        assert str(e_info.value) == "The specified group-by categories (['invalid']) were not found."

    def test_filter_groupby_invalid_warn(self, valid_metadata: pd.DataFrame, capsys):
        groups = ['country', 'year', 'month', 'invalid']
        metadata = valid_metadata.copy()
        strains = metadata.index.tolist()
        group_by_strain = augur.filter.subsample.get_groups_for_subsampling(strains, metadata, group_by=groups)
        assert group_by_strain == {
            'SEQ_1': ('A', 2020, '2020-01', 'unknown'),
            'SEQ_2': ('A', 2020, '2020-02', 'unknown'),
            'SEQ_3': ('B', 2020, '2020-03', 'unknown'),
            'SEQ_4': ('B', 2020, '2020-04', 'unknown'),
            'SEQ_5': ('B', 2020, '2020-05', 'unknown')
        }
        captured = capsys.readouterr()
        assert captured.err == "WARNING: Some of the specified group-by categories couldn't be found: invalid\nFiltering by group may behave differently than expected!\n"

    def test_filter_groupby_missing_year_error(self, valid_metadata: pd.DataFrame):
        groups = ['year']
        metadata = valid_metadata.copy()
        metadata = metadata.drop('date', axis='columns')
        strains = metadata.index.tolist()
        with pytest.raises(AugurError) as e_info:
            augur.filter.subsample.get_groups_for_subsampling(strains, metadata, group_by=groups)
        assert str(e_info.value) == "The specified group-by categories (['year']) were not found. Note that using any of ['month', 'week', 'year'] requires a column called 'date'."

    def test_filter_groupby_missing_month_error(self, valid_metadata: pd.DataFrame):
        groups = ['month']
        metadata = valid_metadata.copy()
        metadata = metadata.drop('date', axis='columns')
        strains = metadata.index.tolist()
        with pytest.raises(AugurError) as e_info:
            augur.filter.subsample.get_groups_for_subsampling(strains, metadata, group_by=groups)
        assert str(e_info.value) == "The specified group-by categories (['month']) were not found. Note that using any of ['month', 'week', 'year'] requires a column called 'date'."

    def test_filter_groupby_missing_year_and_month_error(self, valid_metadata: pd.DataFrame):
        groups = ['year', 'month']
        metadata = valid_metadata.copy()
        metadata = metadata.drop('date', axis='columns')
        strains = metadata.index.tolist()
        with pytest.raises(AugurError) as e_info:
            augur.filter.subsample.get_groups_for_subsampling(strains, metadata, group_by=groups)
        assert str(e_info.value) == "The specified group-by categories (['year', 'month']) were not found. Note that using any of ['month', 'week', 'year'] requires a column called 'date'."

    def test_filter_groupby_missing_date_warn(self, valid_metadata: pd.DataFrame, capsys):
        groups = ['country', 'year', 'month']
        metadata = valid_metadata.copy()
        metadata = metadata.drop('date', axis='columns')
        strains = metadata.index.tolist()
        group_by_strain = augur.filter.subsample.get_groups_for_subsampling(strains, metadata, group_by=groups)
        assert group_by_strain == {
            'SEQ_1': ('A', 'unknown', 'unknown'),
            'SEQ_2': ('A', 'unknown', 'unknown'),
            'SEQ_3': ('B', 'unknown', 'unknown'),
            'SEQ_4': ('B', 'unknown', 'unknown'),
            'SEQ_5': ('B', 'unknown', 'unknown')
        }
        captured = capsys.readouterr()
        assert captured.err == "WARNING: A 'date' column could not be found to group-by ['month', 'year'].\nFiltering by group may behave differently than expected!\n"

    def test_filter_groupby_no_strains(self, valid_metadata: pd.DataFrame):
        groups = ['country', 'year', 'month']
        metadata = valid_metadata.copy()
        strains = []
        group_by_strain = augur.filter.subsample.get_groups_for_subsampling(strains, metadata, group_by=groups)
        assert group_by_strain == {}

    def test_filter_groupby_only_year_provided(self, valid_metadata: pd.DataFrame):
        groups = ['country', 'year']
        metadata = valid_metadata.copy()
        metadata['date'] = '2020'
        strains = metadata.index.tolist()
        group_by_strain = augur.filter.subsample.get_groups_for_subsampling(strains, metadata, group_by=groups)
        assert group_by_strain == {
            'SEQ_1': ('A', 2020),
            'SEQ_2': ('A', 2020),
            'SEQ_3': ('B', 2020),
            'SEQ_4': ('B', 2020),
            'SEQ_5': ('B', 2020)
        }

    def test_filter_groupby_only_year_month_provided(self, valid_metadata: pd.DataFrame):
        groups = ['country', 'year', 'month']
        metadata = valid_metadata.copy()
        metadata['date'] = '2020-01'
        strains = metadata.index.tolist()
        group_by_strain = augur.filter.subsample.get_groups_for_subsampling(strains, metadata, group_by=groups)
        assert group_by_strain == {
            'SEQ_1': ('A', 2020, '2020-01'),
            'SEQ_2': ('A', 2020, '2020-01'),
            'SEQ_3': ('B', 2020, '2020-01'),
            'SEQ_4': ('B', 2020, '2020-01'),
            'SEQ_5': ('B', 2020, '2020-01')
        }

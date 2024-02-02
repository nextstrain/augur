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
    def test_filter_groupby_invalid_error(self):
        groups = ['invalid']
        metadata_columns = {'strain', 'date', 'country'}
        with pytest.raises(AugurError) as e_info:
            augur.filter.subsample.get_valid_group_by_columns(metadata_columns, groups)
        assert str(e_info.value) == "The specified group-by categories (['invalid']) were not found."

    def test_filter_groupby_invalid_warn(self, capsys):
        groups = ['country', 'year', 'month', 'invalid']
        metadata_columns = {'strain', 'date', 'country'}
        augur.filter.subsample.get_valid_group_by_columns(metadata_columns, groups)
        captured = capsys.readouterr()
        assert captured.err == "WARNING: Some of the specified group-by categories couldn't be found: invalid\nFiltering by group may behave differently than expected!\n"

    def test_filter_groupby_missing_year_error(self):
        groups = ['year']
        metadata_columns = {'strain', 'country'}
        with pytest.raises(AugurError) as e_info:
            augur.filter.subsample.get_valid_group_by_columns(metadata_columns, groups)
        assert str(e_info.value) == "The specified group-by categories (['year']) were not found. Note that using any of ['month', 'week', 'year'] requires a column called 'date'."

    def test_filter_groupby_missing_month_error(self):
        groups = ['month']
        metadata_columns = {'strain', 'country'}
        with pytest.raises(AugurError) as e_info:
            augur.filter.subsample.get_valid_group_by_columns(metadata_columns, groups)
        assert str(e_info.value) == "The specified group-by categories (['month']) were not found. Note that using any of ['month', 'week', 'year'] requires a column called 'date'."

    def test_filter_groupby_missing_year_and_month_error(self):
        groups = ['year', 'month']
        metadata_columns = {'strain', 'country'}
        with pytest.raises(AugurError) as e_info:
            augur.filter.subsample.get_valid_group_by_columns(metadata_columns, groups)
        assert str(e_info.value) == "The specified group-by categories (['year', 'month']) were not found. Note that using any of ['month', 'week', 'year'] requires a column called 'date'."

    def test_filter_groupby_missing_date_warn(self, capsys):
        groups = ['country', 'year', 'month']
        metadata_columns = {'strain', 'country'}
        augur.filter.subsample.get_valid_group_by_columns(metadata_columns, groups)
        captured = capsys.readouterr()
        assert captured.err == "WARNING: A 'date' column could not be found to group-by ['month', 'year'].\nFiltering by group may behave differently than expected!\n"

import pytest
import pandas as pd
from augur.filter import get_groups_for_subsampling, FilterException

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

class TestFilterGroupBy:
    def test_filter_groupby_strain_subset(self, valid_metadata: pd.DataFrame):
        metadata = valid_metadata.copy()
        strains = ['SEQ_1', 'SEQ_3', 'SEQ_5']
        group_by_strain, skipped_strains = get_groups_for_subsampling(strains, metadata)
        assert group_by_strain == {
            'SEQ_1': ('_dummy',),
            'SEQ_3': ('_dummy',),
            'SEQ_5': ('_dummy',)
        }
        assert skipped_strains == []

    def test_filter_groupby_dummy(self, valid_metadata: pd.DataFrame):
        metadata = valid_metadata.copy()
        strains = metadata.index.tolist()
        group_by_strain, skipped_strains = get_groups_for_subsampling(strains, metadata)
        assert group_by_strain == {
            'SEQ_1': ('_dummy',),
            'SEQ_2': ('_dummy',),
            'SEQ_3': ('_dummy',),
            'SEQ_4': ('_dummy',),
            'SEQ_5': ('_dummy',)
        }
        assert skipped_strains == []

    def test_filter_groupby_invalid_error(self, valid_metadata: pd.DataFrame):
        groups = ['invalid']
        metadata = valid_metadata.copy()
        strains = metadata.index.tolist()
        with pytest.raises(FilterException) as e_info:
            get_groups_for_subsampling(strains, metadata, group_by=groups)
        assert str(e_info.value) == "The specified group-by categories (['invalid']) were not found. No sequences-per-group sampling will be done."

    def test_filter_groupby_invalid_warn(self, valid_metadata: pd.DataFrame, capsys):
        groups = ['country', 'year', 'month', 'invalid']
        metadata = valid_metadata.copy()
        strains = metadata.index.tolist()
        group_by_strain, _ = get_groups_for_subsampling(strains, metadata, group_by=groups)
        assert group_by_strain == {
            'SEQ_1': ('A', 2020, (2020, 1), 'unknown'),
            'SEQ_2': ('A', 2020, (2020, 2), 'unknown'),
            'SEQ_3': ('B', 2020, (2020, 3), 'unknown'),
            'SEQ_4': ('B', 2020, (2020, 4), 'unknown'),
            'SEQ_5': ('B', 2020, (2020, 5), 'unknown')
        }
        captured = capsys.readouterr()
        assert captured.err == "WARNING: Some of the specified group-by categories couldn't be found: invalid\nFiltering by group may behave differently than expected!\n"

    def test_filter_groupby_skip_year(self, valid_metadata: pd.DataFrame):
        groups = ['country', 'year', 'month']
        metadata = valid_metadata.copy()
        metadata.at["SEQ_2", "date"] = "XXXX-02-01"
        strains = metadata.index.tolist()
        group_by_strain, skipped_strains = get_groups_for_subsampling(strains, metadata, group_by=groups)
        assert group_by_strain == {
            'SEQ_1': ('A', 2020, (2020, 1)),
            'SEQ_3': ('B', 2020, (2020, 3)),
            'SEQ_4': ('B', 2020, (2020, 4)),
            'SEQ_5': ('B', 2020, (2020, 5))
        }
        assert skipped_strains == [{'strain': 'SEQ_2', 'filter': 'skip_group_by_with_ambiguous_year', 'kwargs': ''}]

    def test_filter_groupby_skip_month(self, valid_metadata: pd.DataFrame):
        groups = ['country', 'year', 'month']
        metadata = valid_metadata.copy()
        metadata.at["SEQ_2", "date"] = "2020-XX-01"
        strains = metadata.index.tolist()
        group_by_strain, skipped_strains = get_groups_for_subsampling(strains, metadata, group_by=groups)
        assert group_by_strain == {
            'SEQ_1': ('A', 2020, (2020, 1)),
            'SEQ_3': ('B', 2020, (2020, 3)),
            'SEQ_4': ('B', 2020, (2020, 4)),
            'SEQ_5': ('B', 2020, (2020, 5))
        }
        assert skipped_strains == [{'strain': 'SEQ_2', 'filter': 'skip_group_by_with_ambiguous_month', 'kwargs': ''}]

    def test_filter_groupby_skip_month_2(self, valid_metadata: pd.DataFrame):
        groups = ['country', 'year', 'month']
        metadata = valid_metadata.copy()
        metadata.at["SEQ_2", "date"] = "2020"
        strains = metadata.index.tolist()
        group_by_strain, skipped_strains = get_groups_for_subsampling(strains, metadata, group_by=groups)
        assert group_by_strain == {
            'SEQ_1': ('A', 2020, (2020, 1)),
            'SEQ_3': ('B', 2020, (2020, 3)),
            'SEQ_4': ('B', 2020, (2020, 4)),
            'SEQ_5': ('B', 2020, (2020, 5))
        }
        assert skipped_strains == [{'strain': 'SEQ_2', 'filter': 'skip_group_by_with_ambiguous_month', 'kwargs': ''}]

    def test_filter_groupby_missing_year_error(self, valid_metadata: pd.DataFrame):
        groups = ['year']
        metadata = valid_metadata.copy()
        metadata = metadata.drop('date', axis='columns')
        strains = metadata.index.tolist()
        with pytest.raises(FilterException) as e_info:
            get_groups_for_subsampling(strains, metadata, group_by=groups)
        assert str(e_info.value) == "The specified group-by categories (['year']) were not found. No sequences-per-group sampling will be done. Note that using 'year' or 'year month' requires a column called 'date'."

    def test_filter_groupby_missing_month_error(self, valid_metadata: pd.DataFrame):
        groups = ['month']
        metadata = valid_metadata.copy()
        metadata = metadata.drop('date', axis='columns')
        strains = metadata.index.tolist()
        with pytest.raises(FilterException) as e_info:
            get_groups_for_subsampling(strains, metadata, group_by=groups)
        assert str(e_info.value) == "The specified group-by categories (['month']) were not found. No sequences-per-group sampling will be done. Note that using 'year' or 'year month' requires a column called 'date'."

    def test_filter_groupby_missing_year_and_month_error(self, valid_metadata: pd.DataFrame):
        groups = ['year', 'month']
        metadata = valid_metadata.copy()
        metadata = metadata.drop('date', axis='columns')
        strains = metadata.index.tolist()
        with pytest.raises(FilterException) as e_info:
            get_groups_for_subsampling(strains, metadata, group_by=groups)
        assert str(e_info.value) == "The specified group-by categories (['year', 'month']) were not found. No sequences-per-group sampling will be done. Note that using 'year' or 'year month' requires a column called 'date'."

    def test_filter_groupby_missing_date_warn(self, valid_metadata: pd.DataFrame, capsys):
        groups = ['country', 'year', 'month']
        metadata = valid_metadata.copy()
        metadata = metadata.drop('date', axis='columns')
        strains = metadata.index.tolist()
        group_by_strain, skipped_strains = get_groups_for_subsampling(strains, metadata, group_by=groups)
        assert group_by_strain == {
            'SEQ_1': ('A', 'unknown', 'unknown'),
            'SEQ_2': ('A', 'unknown', 'unknown'),
            'SEQ_3': ('B', 'unknown', 'unknown'),
            'SEQ_4': ('B', 'unknown', 'unknown'),
            'SEQ_5': ('B', 'unknown', 'unknown')
        }
        captured = capsys.readouterr()
        assert captured.err == "WARNING: A 'date' column could not be found to group-by year or month.\nFiltering by group may behave differently than expected!\n"
        assert skipped_strains == []

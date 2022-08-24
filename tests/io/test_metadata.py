import pytest

from augur.errors import AugurError
from augur.io.metadata import read_table_to_dict


@pytest.fixture
def expected_record():
    return {
        'strain': 'SEQ_A',
        'date': '2020-10-03',
        'country': 'USA'
    }

class TestReadMetadataToDict:
    def test_read_table_to_dict_with_csv(self, tmpdir, expected_record):
        path = str(tmpdir / 'metadata.csv')
        with open(path, 'w') as fh:
            fh.write('strain,date,country\n')
            fh.write('SEQ_A,2020-10-03,USA\n')

        record = next(read_table_to_dict(path))
        assert record == expected_record

    def test_read_table_to_dict_with_tsv(self, tmpdir, expected_record):
        path = str(tmpdir / 'metadata.tsv')
        with open(path, 'w') as fh:
            fh.write('strain\tdate\tcountry\n')
            fh.write('SEQ_A\t2020-10-03\tUSA\n')

        record = next(read_table_to_dict(path))
        assert record == expected_record

    def test_read_table_to_dict_with_bad_delimiter(self, tmpdir):
        path = str(tmpdir / 'metadata.txt')
        with open(path, 'w') as fh:
            fh.write('strain date country\n')
            fh.write('SEQ_A 2020-10-03 USA\n')

        with pytest.raises(AugurError) as e_info:
            next(read_table_to_dict(path))

        assert str(e_info.value) == f"Could not determine the delimiter of {path!r}. File must be a CSV or TSV."

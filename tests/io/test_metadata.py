import pytest

from augur.errors import AugurError
from augur.io.metadata import read_table_to_dict
from augur.types import DataErrorMethod


@pytest.fixture
def expected_record():
    return {
        'strain': 'SEQ_A',
        'date': '2020-10-03',
        'country': 'USA'
    }

@pytest.fixture
def metadata_with_duplicate(tmpdir):
    path = str(tmpdir / 'metadata.tsv')
    with open(path, 'w') as fh:
        fh.write('strain\tdate\tcountry\n')
        fh.write('SEQ_A\t2020-10-03\tUSA\n')
        fh.write('SEQ_A\t2020-10-03\tUSA\n')
        fh.write('SEQ_B\t2020-10-03\tUSA\n')
        fh.write('SEQ_B\t2020-10-03\tUSA\n')
    return path

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

    @pytest.mark.parametrize('id_column', ['strain', None])
    def test_read_table_to_dict_with_duplicates(self, metadata_with_duplicate, id_column):
        with pytest.raises(AugurError) as e_info:
            list(read_table_to_dict(metadata_with_duplicate, id_column=id_column))
        assert str(e_info.value) == f"Encountered record with duplicate id 'SEQ_A' in {metadata_with_duplicate!r}"

    @pytest.mark.parametrize('id_column', ['strain', None])
    def test_read_table_to_dict_with_duplicates_error_all(self, metadata_with_duplicate, id_column):
        with pytest.raises(AugurError) as e_info:
            list(read_table_to_dict(metadata_with_duplicate, DataErrorMethod("error_all"), id_column=id_column))
        assert str(e_info.value) == f"The following records are duplicated in {metadata_with_duplicate!r}:\n'SEQ_A'\n'SEQ_B'"

    @pytest.mark.parametrize('id_column', ['strain', None])
    def test_read_table_to_dict_with_duplicates_warning(self, capsys, metadata_with_duplicate, id_column):
        list(read_table_to_dict(metadata_with_duplicate, DataErrorMethod('warn'), id_column=id_column))
        captured = capsys.readouterr()
        assert captured.err == (
            f"WARNING: Encountered record with duplicate id 'SEQ_A' in {metadata_with_duplicate!r}\n"
            f"WARNING: Encountered record with duplicate id 'SEQ_B' in {metadata_with_duplicate!r}\n"
            f"WARNING: The following records are duplicated in {metadata_with_duplicate!r}:\n'SEQ_A'\n'SEQ_B'\n"
        )

    def test_read_table_to_dict_with_duplicates_silent(self, capsys, metadata_with_duplicate):
        list(read_table_to_dict(metadata_with_duplicate, DataErrorMethod('silent')))
        assert "WARNING" not in capsys.readouterr().err

    def test_read_table_to_dict_with_duplicate_and_bad_id(self, metadata_with_duplicate):
        id_column = "bad_id"
        with pytest.raises(AugurError) as e_info:
            list(read_table_to_dict(metadata_with_duplicate, id_column=id_column))
        assert str(e_info.value) == f"The provided id column {id_column!r} does not exist in {metadata_with_duplicate!r}."

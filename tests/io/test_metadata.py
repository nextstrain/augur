import pytest
import shutil

from augur.errors import AugurError
from augur.io.metadata import read_table_to_dict, read_metadata_with_sequences
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


@pytest.fixture
def fasta_file(tmpdir):
    path = str(tmpdir / 'sequences.fasta')
    with open(path, 'w') as fh:
        fh.writelines([
            '>SEQ_A\nAAAA\n',
            '>SEQ_T\nTTTT\n',
            '>SEQ_C\nCCCC\n',
            '>SEQ_G\nGGGG\n'
        ])
    return path

@pytest.fixture
def metadata_file(tmpdir):
    path = str(tmpdir / 'metadata.tsv')
    with open(path, 'w') as fh:
        fh.writelines([
            'strain\tcountry\tdate\n',
            'SEQ_A\tUSA\t2020-10-01\n',
            'SEQ_T\tUSA\t2020-10-02\n',
            'SEQ_C\tUSA\t2020-10-03\n',
            'SEQ_G\tUSA\t2020-10-04\n'
        ])
    return path

def unmatched_sequences():
    return [
        '>EXTRA_SEQ_A\nAAAAA\n',
        '>EXTRA_SEQ_T\nTTTTT\n'
    ]

def unmatched_metadata():
    return [
        'EXTRA_METADATA_A\tUSA\t2020-10-01\n',
        'EXTRA_METADATA_T\tUSA\t2020-10-02\n',
    ]

def copy_and_append_to_file(src, dst, appended_content):
    shutil.copy(src, dst)
    with open(dst, 'a') as fh:
        fh.writelines(appended_content)
    return dst

@pytest.fixture
def fasta_with_unmatched(tmpdir, fasta_file):
    path = str(tmpdir / 'extra-sequences.fasta')
    return copy_and_append_to_file(fasta_file, path, unmatched_sequences())

@pytest.fixture
def metadata_with_unmatched(tmpdir, metadata_file):
    path = str(tmpdir / 'extra-metadata.tsv')
    return copy_and_append_to_file(metadata_file, path, unmatched_metadata())

class TestReadMetadataWithSequence:
    def test_read_metadata_with_sequence(self, metadata_file, fasta_file):
        records = list(read_metadata_with_sequences(metadata_file, fasta_file, 'strain'))
        assert len(records) == 4
        for record in records:
            seq_base = record['strain'].split("_")[-1].upper()
            expected_sequence = seq_base * 4
            assert record['sequence'] == expected_sequence

    def test_read_metadata_with_sequences_with_bad_id(self, metadata_file, fasta_file):
        id_field = "bad_id"
        with pytest.raises(AugurError) as e_info:
            next(read_metadata_with_sequences(metadata_file, fasta_file, id_field))
        assert str(e_info.value) == f"The provided sequence id column {id_field!r} does not exist in the metadata."

    def test_read_metadata_with_sequences_with_unmatched(self, metadata_with_unmatched, fasta_with_unmatched):
        with pytest.raises(AugurError) as e_info:
            list(read_metadata_with_sequences(metadata_with_unmatched, fasta_with_unmatched, 'strain'))
        assert str(e_info.value) == "Encountered metadata record 'EXTRA_METADATA_A' without a matching sequence."

    def test_read_metadata_with_sequences_with_unmatched_error_all(self, metadata_with_unmatched, fasta_with_unmatched):
        with pytest.raises(AugurError) as e_info:
            list(read_metadata_with_sequences(
                metadata_with_unmatched,
                fasta_with_unmatched,
                'strain',
                unmatched_reporting=DataErrorMethod.ERROR_ALL))
        assert str(e_info.value) == (
            "The output may be incomplete because there were unmatched records.\n"
            "The following metadata records did not have a matching sequence:\n"
            "'EXTRA_METADATA_A'\n'EXTRA_METADATA_T'\n"
            "The following sequence records did not have a matching metadata record:\n"
            "'EXTRA_SEQ_A'\n'EXTRA_SEQ_T'"
        )

    def test_read_metadata_with_sequences_with_unmatched_warning(self, capsys, metadata_with_unmatched, fasta_with_unmatched):
        records = list(read_metadata_with_sequences(
            metadata_with_unmatched,
            fasta_with_unmatched,
            'strain',
            unmatched_reporting=DataErrorMethod.WARN))
        assert len(records) == 4
        assert [record['strain'] for record in records] == ['SEQ_A', 'SEQ_T', 'SEQ_C', 'SEQ_G']

        captured = capsys.readouterr()
        assert captured.err == (
            "WARNING: Encountered metadata record 'EXTRA_METADATA_A' without a matching sequence.\n"
            "WARNING: Encountered metadata record 'EXTRA_METADATA_T' without a matching sequence.\n"
            "WARNING: The output may be incomplete because there were unmatched records.\n"
            "The following metadata records did not have a matching sequence:\n"
            "'EXTRA_METADATA_A'\n'EXTRA_METADATA_T'\n"
            "The following sequence records did not have a matching metadata record:\n"
            "'EXTRA_SEQ_A'\n'EXTRA_SEQ_T'\n"
        )

    def test_read_metadata_with_sequences_with_unmatched_silent(self, capsys, metadata_with_unmatched, fasta_with_unmatched):
        records = list(read_metadata_with_sequences(
            metadata_with_unmatched,
            fasta_with_unmatched,
            'strain',
            unmatched_reporting=DataErrorMethod.SILENT))
        assert len(records) == 4
        assert [record['strain'] for record in records] == ['SEQ_A', 'SEQ_T', 'SEQ_C', 'SEQ_G']
        assert "WARNING" not in capsys.readouterr().err

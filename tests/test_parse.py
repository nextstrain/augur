import Bio.SeqIO
import Bio.Seq
from pathlib import Path
import sys

# make sure we get the local version of parse
sys.path.insert(0, (str(Path(__file__).parent.parent.parent)))
from augur import parse


class TestParse:
    def test_fix_dates(self, capsys):
        full_date = "4-5-2020"
        assert parse.fix_dates(full_date) == "2020-05-04"
        assert parse.fix_dates(full_date, dayfirst=False) == "2020-04-05"
        partial_date_no_month = "2020-04"
        assert parse.fix_dates(partial_date_no_month) == "2020-04-XX"
        assert parse.fix_dates(partial_date_no_month, dayfirst=False) == "2020-04-XX"

        partial_date_only_year = "2020"
        assert parse.fix_dates(partial_date_only_year) == "2020-XX-XX"

        # We should allow valid ambiguous dates to pass through the function
        # without a warning to the user.
        valid_ambiguous_date = "2020-XX-XX"
        assert parse.fix_dates(valid_ambiguous_date) == valid_ambiguous_date
        assert "WARNING" not in capsys.readouterr().err

        # We should warn the user when their dates cannot be parsed by pandas or
        # as a valid ambiguous date.
        malformed_date = "Aaasd123AS"
        assert(parse.fix_dates(malformed_date)) == malformed_date
        assert "WARNING" in capsys.readouterr().err

    def test_prettify(self):
        # test trim
        trim_string = "trim_here"
        assert parse.prettify(trim_string, trim=4) == trim_string[:4] + "..."
        trim_string_short = "trim"
        assert parse.prettify(trim_string_short, trim=4) == trim_string_short

        # test caps
        caps_string_pos = "usa"
        caps_string_neg = "usad"
        assert parse.prettify(caps_string_pos) == caps_string_pos.upper()
        assert parse.prettify(caps_string_neg) == caps_string_neg

        # test camelCase
        camel_case_string = "testing_one_two_three"
        assert parse.prettify(camel_case_string, camelCase=True) == "Testing One Two Three"

        # test remove comma
        remove_comma_string = "remove,the,commas"
        parse.prettify(remove_comma_string, removeComma=True)
        assert parse.prettify(remove_comma_string, removeComma=True) == "removethecommas"

        # test etal
        etal_lower_string = "testing string Et Al Et al"
        etal_strip_string = "nextstrain et al. et al Et Al."
        assert parse.prettify(etal_lower_string, etal='lower') == etal_lower_string.lower()
        assert parse.prettify(etal_strip_string, etal='strip') == "nextstrain   "

    def test_parse_sequence(self):
        sequence_record = Bio.SeqRecord.SeqRecord(
            seq=Bio.Seq.Seq("ATCG"),
            id="sequence_A",
            description="sequence_A|USA|2020-10-03|north_america"
        )
        fields = ["strain", "country", "date", "region"]
        sequence, metadata = parse.parse_sequence(
            sequence_record,
            fields=fields,
            strain_key="strain",
            separator="|",
            prettify_fields=["region"],
            fix_dates_format=None,
        )

        assert sequence.id == metadata["strain"]
        assert sequence.id == "sequence_A"
        assert set(fields) == set(metadata.keys())
        assert metadata["region"] == "North America"

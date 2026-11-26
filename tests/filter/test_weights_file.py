import pytest

from augur.filter.subsample import get_weighted_group_sizes
from augur.filter.weights_file import InvalidWeightsFile, get_default_weight, read_weights_file


class TestReadWeightsFile:
    @pytest.mark.parametrize("literal_group", ["NA", "N/A", "None", "NULL", "nan", ""])
    def test_preserves_literal_na_tokens_as_strings(self, tmp_path, literal_group):
        """Match read_metadata(na_filter=False) so tokens like 'NA' (Namibia / influenza NA) stay literal strings."""
        weights_file = tmp_path / "weights.tsv"
        weights_file.write_text(f"country_code\tweight\n{literal_group}\t3\nZA\t1\n", encoding="utf-8")

        weights = read_weights_file(str(weights_file))

        assert weights["country_code"].tolist() == [literal_group, "ZA"]
        assert weights["weight"].tolist() == [3, 1]

    @pytest.mark.parametrize("bad_weight", ["NA", "None", ""])
    def test_rejects_na_or_empty_weight_values(self, tmp_path, bad_weight):
        weights_file = tmp_path / "weights.tsv"
        weights_file.write_text(f"country_code\tweight\nZA\t{bad_weight}\n", encoding="utf-8")

        with pytest.raises(InvalidWeightsFile, match=r"Found non-numeric weights on the following lines: \[2\]"):
            read_weights_file(str(weights_file))

    def test_partial_default_row_with_empty_column_is_not_default(self, tmp_path):
        weights_file = tmp_path / "weights.tsv"
        weights_file.write_text(
            "location\tmonth\tweight\nA\t2000-01\t2\ndefault\t\t1\n",
            encoding="utf-8",
        )

        weights = read_weights_file(str(weights_file))

        assert get_default_weight(weights, ["location", "month"]) is None


class TestWeightedGroupSizesLiteralNA:
    def test_explicit_weight_for_literal_na_not_overridden_by_default(self, tmp_path):
        weights_file = tmp_path / "weights.tsv"
        weights_file.write_text(
            "country_code\tweight\nNA\t3\nZA\t1\ndefault\t1\n",
            encoding="utf-8",
        )

        sizes = get_weighted_group_sizes(
            records_per_group={("NA",): 4, ("ZA",): 4},
            group_by=["country_code"],
            weights_file=str(weights_file),
            target_total_size=4,
            output_sizes_file=None,
            random_seed=42,
        )

        assert sizes == {("NA",): 3, ("ZA",): 1}

    def test_explicit_zero_weight_excludes_literal_na_when_default_present(self, tmp_path):
        weights_file = tmp_path / "weights.tsv"
        weights_file.write_text(
            "country_code\tweight\nNA\t0\nZA\t1\ndefault\t1\n",
            encoding="utf-8",
        )

        sizes = get_weighted_group_sizes(
            records_per_group={("NA",): 4, ("ZA",): 4},
            group_by=["country_code"],
            weights_file=str(weights_file),
            target_total_size=4,
            output_sizes_file=None,
            random_seed=42,
        )

        assert sizes == {("NA",): 0, ("ZA",): 4}

    def test_complete_weights_with_literal_na_without_default(self, tmp_path):
        weights_file = tmp_path / "weights.tsv"
        weights_file.write_text(
            "country_code\tweight\nNA\t3\nZA\t1\n",
            encoding="utf-8",
        )

        sizes = get_weighted_group_sizes(
            records_per_group={("NA",): 4, ("ZA",): 4},
            group_by=["country_code"],
            weights_file=str(weights_file),
            target_total_size=4,
            output_sizes_file=None,
            random_seed=42,
        )

        assert sizes == {("NA",): 3, ("ZA",): 1}

    def test_literal_na_with_unweighted_column(self, tmp_path):
        weights_file = tmp_path / "weights.tsv"
        weights_file.write_text(
            "country_code\tweight\nNA\t3\nZA\t1\ndefault\t1\n",
            encoding="utf-8",
        )

        sizes = get_weighted_group_sizes(
            records_per_group={
                ("NA", "2024"): 4,
                ("NA", "2025"): 4,
                ("ZA", "2024"): 4,
                ("ZA", "2025"): 4,
            },
            group_by=["country_code", "year"],
            weights_file=str(weights_file),
            target_total_size=8,
            output_sizes_file=None,
            random_seed=42,
        )

        assert sizes == {
            ("NA", "2024"): 3,
            ("NA", "2025"): 3,
            ("ZA", "2024"): 1,
            ("ZA", "2025"): 1,
        }

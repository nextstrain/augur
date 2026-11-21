import argparse
from pathlib import Path
import pytest

from augur.filter import register_arguments
from augur.filter._run import run


def _write_inputs(
    tmp_path: Path,
    metadata_strains: list[str],
    index_strains: list[str],
) -> tuple[Path, Path, Path, Path]:
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "strain\tdate\n" + "".join(f"{s}\t2024-01-01\n" for s in metadata_strains),
        encoding="utf-8",
    )

    sequences = tmp_path / "sequences.fasta"
    sequences.write_text(">placeholder\nATCG\n", encoding="utf-8")

    sequence_index = tmp_path / "sequence_index.tsv"
    sequence_index.write_text(
        "strain\tlength\tA\tC\tG\tT\tN\tother_IUPAC\t?\t-\tinvalid_nucleotides\n"
        + "".join(f"{s}\t4\t1\t1\t1\t1\t0\t0\t0\t0\t0\n" for s in index_strains),
        encoding="utf-8",
    )

    output_strains = tmp_path / "filtered_strains.txt"
    return metadata, sequences, sequence_index, output_strains


@pytest.mark.parametrize(
    "sequence_strains, index_strains, expect_warning",
    [
        # Partial overlap ({A, B} vs {A, C}): previously missed by `>` check (#1825)
        (["A", "B"], ["A", "C"], True),
        # Strict superset ({A, B} vs {A}): must warn
        (["A", "B"], ["A"], True),
        # Disjoint ({A, B} vs {C}): must warn
        (["A", "B"], ["C"], True),
        # Equal ({A, B} vs {A, B}): no warning
        (["A", "B"], ["A", "B"], False),
        # Proper subset ({A} vs {A, B}): sequence index is a superset, no warning
        (["A"], ["A", "B"], False),
    ],
)
def test_sequence_index_mismatch_warning(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    sequence_strains: list[str],
    index_strains: list[str],
    expect_warning: bool,
) -> None:
    metadata, sequences, sequence_index, output_strains = _write_inputs(
        tmp_path,
        metadata_strains=["A", "B", "C"],
        index_strains=index_strains,
    )
    monkeypatch.setattr(
        "augur.filter._run.read_sequence_ids",
        lambda _path, _nthreads: set(sequence_strains),
    )

    parser = argparse.ArgumentParser()
    register_arguments(parser)
    args = parser.parse_args(
        [
            "--metadata",
            str(metadata),
            "--sequences",
            str(sequences),
            "--sequence-index",
            str(sequence_index),
            "--output-strains",
            str(output_strains),
        ]
    )

    run(args)

    captured = capsys.readouterr()
    warning = "WARNING: The sequence index is out of sync with the provided sequences."
    if expect_warning:
        assert warning in captured.err
    else:
        assert warning not in captured.err

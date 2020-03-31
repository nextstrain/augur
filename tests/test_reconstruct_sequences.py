import pytest
from Bio.Align import MultipleSeqAlignment
from augur.reconstruct_sequences import (
    get_sequence,
    load_alignments
)


PARENT_SEQUENCE = "SAWT"


@pytest.mark.parametrize(
    'parent,muts,output',
    (
        (PARENT_SEQUENCE, ['T4S'], 'SAWS'),          # Single mutation
        (PARENT_SEQUENCE, ['S1W', 'T4S'], 'WAWS')    # Double mutation
    )
)
def test_get_sequence(parent, muts, output):
    out = get_sequence(parent, muts)
    assert out == output


@pytest.mark.parametrize(
    'parent,muts',
    (
        (PARENT_SEQUENCE, ['T5S']),          # Bad index at site 5
        (PARENT_SEQUENCE, ['S1W', 'T5S'])    # Double mutation, bad index at site 5
    )
)
def test_get_sequences_index_error(parent, muts):
    with pytest.raises(IndexError) as error_info:
        get_sequence(parent, muts)


@pytest.mark.parametrize(
    'parent,muts',
    (
        (PARENT_SEQUENCE, ['A4S']),          # Bad index at site 5
        (PARENT_SEQUENCE, ['S1W', 'A4S'])    # Double mutation, bad index at site 5
    )
)
def test_get_sequences_assert_error(parent, muts):
    with pytest.raises(AssertionError) as error_info:
        get_sequence(parent, muts)

    assert "does not match" in str(error_info.value)


ALIGNMENT1 = """\
>seq1
AKWTKLMNQPKA-
>seq2
MKTALAMPWQGTM
"""
ALIGNMENT2 = """\
>seq1
WLKTMTLKPKDGM
>seq2
KLWTYMAPYKLMS
"""

@pytest.fixture
def alignment_data(tmp_path):
    """Create some mock alignments."""
    align1 = tmp_path.joinpath("alignment1.fasta")
    align2 = tmp_path.joinpath("alignment2.fasta")
    align1.write_text(ALIGNMENT1)
    align2.write_text(ALIGNMENT2)


def test_load_alignments(tmp_path, alignment_data):
    align1 = tmp_path.joinpath("alignment1.fasta")
    align2 = tmp_path.joinpath("alignment2.fasta")
    sequence_files = [str(align1), str(align2)]
    gene_names = ['geneA', 'geneB']

    alignments = load_alignments(sequence_files, gene_names)

    assert isinstance(alignments, dict)
    assert 'geneA' in alignments
    assert 'geneB' in alignments
    assert isinstance(alignments['geneA'], MultipleSeqAlignment)

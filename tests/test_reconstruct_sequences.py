import pytest
from augur.reconstruct_sequences import get_sequence


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

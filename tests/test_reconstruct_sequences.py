import pytest
from Bio import  Phylo
from Bio.Align import MultipleSeqAlignment
from augur.utils import read_node_data
from augur.reconstruct_sequences import (
    get_sequence,
    load_alignments,
    reconstruct_sequences_from_tree_of_mutations,
    run
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


@pytest.fixture
def zika_tree_file():
    return 'tests/data/zika_example/results/tree.nwk'


@pytest.fixture
def zika_tree(zika_tree_file):
    return Phylo.read(
        zika_tree_file,
        'newick'
    )


@pytest.fixture
def zika_mutations():
    return 'tests/data/zika_example/results/aa_muts.json'


@pytest.fixture
def zika_node_data(zika_tree_file, zika_mutations):
    return read_node_data(zika_mutations, zika_tree_file)


@pytest.fixture
def zika_root_node(zika_tree):
    return zika_tree.root.name


@pytest.fixture
def zika_gene():
    return "ENV"


@pytest.fixture
def args(tmp_path):
    class Args(): pass
    args = Args()
    args.tree = 'tests/data/zika_example/results/tree.nwk'
    args.gene = 'ENV'
    args.mutations = 'tests/data/zika_example/results/aa_muts.json'
    args.vcf_aa_reference = None
    args.internal_nodes = False
    args.output = str(tmp_path.joinpath('reconstructed_sequences.fasta'))
    return args


def test_reconstruct_sequences_from_tree_of_mutations(
        zika_tree,
        zika_node_data,
        zika_root_node,
        zika_gene
    ):
    sequences, is_terminal = reconstruct_sequences_from_tree_of_mutations(
        zika_tree,
        zika_node_data,
        zika_root_node,
        zika_gene
    )
    assert isinstance(sequences, dict)
    assert isinstance(is_terminal, dict)


def test_run(tmp_path, args):
    run(args)
    output = tmp_path.joinpath('reconstructed_sequences.fasta')
    assert output.exists()
    results = output.read_text()
    assert len(results) > 0
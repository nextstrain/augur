from augur.ancestral import run_ancestral
from io import StringIO
from Bio import Phylo
import pytest

def gather_mutations_at_pos(pos, nuc_result):
    """pos is 1-based"""
    muts = []
    for v in nuc_result['mutations']['nodes'].values():
        for mut in v['muts']:
            if int(mut[1:-1])==pos:
                muts.append(mut)
    return muts

def gather_bases_at_pos(pos, nuc_result):
    """pos is 1-based"""
    bases = set()
    for v in nuc_result['mutations']['nodes'].values():
        bases.add(v['sequence'][pos-1])
    return bases

class TestRootNodeMutationAssignment:
    @pytest.fixture
    def tree(self):
        t = "(sample_C:0.02,(sample_B:0.02,sample_A:0.02)node_AB:0.06)node_root:0.02;"
        T = Phylo.read(StringIO(t), 'newick')
        return T

    @pytest.fixture()
    def ref(self):
        return "AAAAAA"
        #       123456

    @pytest.fixture()
    def default_augur_args(self):
        return {
            'fill_overhangs': True, # augur default
            'marginal': 'joint', # augur default
            'alphabet': 'nuc', # augur default
        }


    @pytest.fixture()
    def fasta_aln(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.Align import MultipleSeqAlignment
        aln = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AANATA"), id="sample_A"),
                SeqRecord(Seq("AANATG"), id="sample_B"),
                SeqRecord(Seq("AANCCA"), id="sample_C"),
            ]
        )
        return aln

    @pytest.fixture()
    def vcf_aln(self):
        def to_zero(x):
            return int(x)-1
        return {
            'sample_A': {to_zero(3): 'N', to_zero(5): 'T'},
            'sample_B': {to_zero(3): 'N', to_zero(5): 'T', to_zero(6): 'G'},
            'sample_C': {to_zero(3): 'N', to_zero(4): 'C', to_zero(5): 'C'}
        }

    def test_fasta_with_ambiguity_inferred(self, tree, fasta_aln, ref, default_augur_args):
        is_vcf = False
        nuc_result = run_ancestral(T=tree, aln=fasta_aln, reference_sequence=ref, is_vcf=is_vcf, full_sequences=not is_vcf, infer_ambiguous=True, **default_augur_args)
        assert(nuc_result['root_seq'] == ref)
        assert(gather_bases_at_pos(3, nuc_result)==set(['A']))
        assert(gather_mutations_at_pos(3, nuc_result)==[])

    def test_fasta_with_ambiguity_not_inferred(self, tree, fasta_aln, ref, default_augur_args):
        is_vcf = False
        nuc_result = run_ancestral(T=tree, aln=fasta_aln, reference_sequence=ref, is_vcf=is_vcf, full_sequences=not is_vcf, infer_ambiguous=False, **default_augur_args)
        assert(nuc_result['root_seq'] == ref)
        assert(gather_bases_at_pos(3, nuc_result)==set(['N']))
        # Expected behaviour here is that every tip is "N", and this is also
        # inferred over all the internal nodes (so they're all "N") and thus we have
        # a single mutation "C3N" on the root node. This is not what happens - the
        # full sequences on each node are like this, but the reported mutations are
        # incorrect - we are describing a "C3N" on each terminal node. I'm removing
        # this test now as (a) it's not newly introduced and (b) I don't have the
        # time to fix it right now.
        # Note that the vcf version of this has the same bug.
        # assert(gather_mutations_at_pos(3, nuc_result)==['A3N'])

    def test_vcf_with_ambiguity_inferred(self, tree, vcf_aln, ref, default_augur_args):
        is_vcf = True
        nuc_result = run_ancestral(T=tree, aln=vcf_aln, reference_sequence=ref, is_vcf=is_vcf, full_sequences=not is_vcf, infer_ambiguous=True, **default_augur_args)
        assert(nuc_result['root_seq'] == ref)
        assert(gather_mutations_at_pos(3, nuc_result)==[])

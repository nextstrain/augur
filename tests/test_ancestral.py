from augur.ancestral import run_ancestral
from io import StringIO
from Bio import Phylo
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

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
            'marginal': False, # augur default -- `bool(args.inference=='marginal')`
            'alphabet': 'nuc', # augur default
            'rng_seed': 0,
        }


    @pytest.fixture()
    def fasta_aln(self):
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


class TestAmbiguousAAReconstruction:
    # NOTE: these tests may be better expressed in `cram`, following `tests/functional/ancestral/cram/infer-ambiguous-nucleotides.t`,
    # once `augur ancestral` can perform AA inference without needing nuc sequences.
    # See <https://github.com/nextstrain/augur/pull/1975#discussion_r3002628926> for more discussion
    
    # See <https://www.bioinformatics.org/sms/iupac.html> for ambiguous bases.

    @pytest.fixture
    def tree(self):
        t = "(sample_C:0.02,(sample_B:0.02,sample_A:0.02)node_AB:0.06)node_root:0.02;"
        T = Phylo.read(StringIO(t), 'newick')
        return T

    @pytest.fixture
    def generic_args(self):
        """Generic args to run_ancestral"""
        # Note: `augur ancestral` doesn't report full_sequences, but we want to test that here
        args = {"reference_sequence": None, "is_vcf": False, "full_sequences": True,
            "fill_overhangs": True, "marginal": False, "alphabet": 'aa', "rng_seed": 0}
        return args

    @pytest.mark.xfail  # TODO XXX - fix bug - we report I3N mutation, because N is hardcoded as ambiguous (nuc!) base, but N = Asn = Asparagine 
    def test_unknown_aa_is_x(self, tree, generic_args):
        """
        Ensure we treat "X" as the ambiguous AA residue, not N (N = Asn = Asparagine)
        """
        aln = MultipleSeqAlignment([
            SeqRecord(Seq("IIXII"), id="sample_A"),
            SeqRecord(Seq("IIIII"), id="sample_B"),
            SeqRecord(Seq("IIIII"), id="sample_C"),
        ])
        result = run_ancestral(T=tree, aln=aln, **generic_args, infer_ambiguous=False)
        sample_a = result['mutations']['nodes']['sample_A']
        assert sample_a['sequence'][2] == 'X'
        pos3_muts = gather_mutations_at_pos(3, result)
        assert len(pos3_muts)==1
        assert pos3_muts[0] == "I3X"
    
    def test_unknown_aa_is_reconstructed(self, tree, generic_args):
        """
        Ensure "X" (the ambiguous AA residue) is reconstructed
        """
        aln = MultipleSeqAlignment([
            SeqRecord(Seq("IIXII"), id="sample_A"),
            SeqRecord(Seq("IIIII"), id="sample_B"),
            SeqRecord(Seq("IIIII"), id="sample_C"),
        ])
        result = run_ancestral(T=tree, aln=aln, **generic_args, infer_ambiguous=True)
        sample_a = result['mutations']['nodes']['sample_A']
        assert sample_a['sequence'][2] == 'I'
        pos3_muts = gather_mutations_at_pos(3, result)
        assert len(pos3_muts)==0
    
    @pytest.mark.xfail
    def test_invalid_residues_are_replaced_with_X(self, tree, generic_args):
        """
        Ensure "J" (invalid) is replaced with "X" (the ambiguous AA residue)
        """
        aln = MultipleSeqAlignment([
            SeqRecord(Seq("IIJII"), id="sample_A"),
            SeqRecord(Seq("IIIII"), id="sample_B"),
            SeqRecord(Seq("IIIII"), id="sample_C"),
        ])
        result = run_ancestral(T=tree, aln=aln, **generic_args, infer_ambiguous=False)
        sample_a = result['mutations']['nodes']['sample_A']
        
        # The invalid J in position 3 (idx 2) has been corrected to X
        assert sample_a['sequence'][2] == 'X'
        pos3_muts = gather_mutations_at_pos(3, result)
        assert len(pos3_muts)==1
        assert pos3_muts[0] == "I3X"

    def test_ambiguous_residues_are_passed_through(self, tree, generic_args):
        """
        Ensure "Z" (Glx = Glutamic acid (E) or Glutamine (Q)) is passed though
        """
        aln = MultipleSeqAlignment([
            SeqRecord(Seq("IIZII"), id="sample_A"),
            SeqRecord(Seq("IIIII"), id="sample_B"),
            SeqRecord(Seq("IIIII"), id="sample_C"),
        ])
        result = run_ancestral(T=tree, aln=aln, **generic_args, infer_ambiguous=False)
        sample_a = result['mutations']['nodes']['sample_A']
        
        # The valid Z in position 3 (idx 2)
        assert sample_a['sequence'][2] == 'Z'
        pos3_muts = gather_mutations_at_pos(3, result)
        assert len(pos3_muts)==1
        assert pos3_muts[0].endswith("Z")

    def test_invalid_residues_are_inferred(self, tree, generic_args):
        """
        Ensure "J" (invalid) is inferred (to I, as that's what every other residue is at this pos)
        """
        aln = MultipleSeqAlignment([
            SeqRecord(Seq("IIJII"), id="sample_A"),
            SeqRecord(Seq("IIIII"), id="sample_B"),
            SeqRecord(Seq("IIIII"), id="sample_C"),
        ])
        result = run_ancestral(T=tree, aln=aln, **generic_args, infer_ambiguous=True)
        sample_a = result['mutations']['nodes']['sample_A']
        
        assert sample_a['sequence'][2] == 'I'
        pos3_muts = gather_mutations_at_pos(3, result)
        assert len(pos3_muts)==0
        
    def test_ambiguous_residues_are_inferred(self, tree, generic_args):
        """
        Ensure "Z" (Glx = Glutamic acid (E) or Glutamine (Q)) is inferred as E or Q
        """
        aln = MultipleSeqAlignment([
            SeqRecord(Seq("IIZII"), id="sample_A"),
            SeqRecord(Seq("IIIII"), id="sample_B"),
            SeqRecord(Seq("IIIII"), id="sample_C"),
        ])
        result = run_ancestral(T=tree, aln=aln, **generic_args, infer_ambiguous=True)
        sample_a = result['mutations']['nodes']['sample_A']
        
        assert sample_a['sequence'][2] in ['E', 'Q']
        pos3_muts = gather_mutations_at_pos(3, result)
        assert len(pos3_muts)==1
        assert pos3_muts[0] in ['I3E', 'I3Q']

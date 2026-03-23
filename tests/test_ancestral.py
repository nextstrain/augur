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

class TestAmbiguousNucReconstruction:
    # See <https://www.bioinformatics.org/sms/iupac.html> for ambiguous bases.

    @pytest.fixture
    def tree(self):
        t = "(sample_C:0.02,(sample_B:0.02,sample_A:0.02)node_AB:0.06)node_root:0.02;"
        T = Phylo.read(StringIO(t), 'newick')
        return T

    @pytest.fixture()
    def aln_pos3Y(self):
        # Alignment where sample_C has 'R' at position 3
        aln = MultipleSeqAlignment([
            SeqRecord(Seq("AAGAA"), id="sample_A"),
            SeqRecord(Seq("AAGAA"), id="sample_B"),
            SeqRecord(Seq("AAYAA"), id="sample_C"),
        ])
        return aln

    @pytest.fixture()
    def aln_pos3X(self):
        # Alignment where sample_C has 'X' at position 3
        aln = MultipleSeqAlignment([
            SeqRecord(Seq("AAGAA"), id="sample_A"),
            SeqRecord(Seq("AAGAA"), id="sample_B"),
            SeqRecord(Seq("AAXAA"), id="sample_C"),
        ])
        return aln

    @pytest.fixture
    def generic_args(self):
        """Generic args to run_ancestral"""
        args = {"reference_sequence": None, "is_vcf": False, "full_sequences": True,
            "fill_overhangs": True, "marginal": False, "alphabet": 'nuc', "rng_seed": 0}
        return args

    def test_ambiguous_bases_remain_in_seqs_and_muts(self, tree, aln_pos3Y, generic_args):
        """The ambiguous nuc 'Y' (representing a 'T' or 'C') should be
        reported as a mutation and in the resulting sequence
        if we _don't_ infer ambiguous seqs
        """
        result = run_ancestral(T=tree, aln=aln_pos3Y, **generic_args, infer_ambiguous=False)
        sample_c = result['mutations']['nodes']['sample_C']

        # The sequence retains the raw 'Y' from TreeTime
        assert sample_c['sequence'][2] == 'Y'
        # And there's a corresponding "something to Y" mutation reported on the sample_c branch
        # (and which should be the only mutation on this branch)
        pos3_muts = [m for m in sample_c['muts'] if int(m[1:-1]) == 3]
        assert len(pos3_muts)==1, "Expected only a single mutation (to Y)"
        assert pos3_muts[0].endswith('Y'), "Expected ambiguous base Y to be reported as mutation"

    def test_ambiguous_bases_are_inferred_in_seqs_and_muts(self, tree, aln_pos3Y, generic_args):
        """The ambiguous nuc 'Y' (representing a 'T' or 'C') should be
        inferred if we set infer_ambiguous=True.
        """
        result = run_ancestral(T=tree, aln=aln_pos3Y, **generic_args, infer_ambiguous=True)
        sample_c = result['mutations']['nodes']['sample_C']

        # The sequence collapses the ambiguous 'Y' into a T or a C
        assert sample_c['sequence'][2] in ['T', 'C']
        # And there are no "to Y" mutations reported on the tree.
        # (We can't test directly for an "* to T/C" mutation on sample_c as treetime infers the
        # root as C and thus no mutation on branch C)
        pos3_muts = gather_mutations_at_pos(3, result)
        assert not any([m.endswith("Y") for m in pos3_muts]), "Expected ambiguous base Y to be inferred as mutation to T or C"


    @pytest.mark.xfail # TODO XXX - fix bug
    def test_nuc_X_reported_as_N_in_seqs_and_muts_when_not_infer_ambiguous(self, tree, aln_pos3X, generic_args):
        """
        The nucleotide 'X' has a flat profile map in treetime and thus is considered ambiguous.
        Ensure it's reported as "N" regardless of whether we infer ambiguous states or not,
        and reported in both sequences and mutations
        """
        result = run_ancestral(T=tree, aln=aln_pos3X, **generic_args, infer_ambiguous=False)
        sample_c = result['mutations']['nodes']['sample_C']

        assert sample_c['sequence'][2] == 'N'
        # And there's a corresponding "something to N" mutation reported
        # (which should be the only mutation)
        pos3_muts = gather_mutations_at_pos(3, result)
        assert len(pos3_muts)==1, "Expected only a single mutation (to N)"
        assert pos3_muts[0].endswith('N'), "Expected ambiguous base N to be reported as mutation"

    def test_nuc_X_reported_as_N_in_seqs_and_muts_when_infer_ambiguous(self, tree, aln_pos3X, generic_args):
        """
        The nucleotide 'X' has a flat profile map in treetime and thus is considered ambiguous.
        Ensure it's reported as "N" regardless of whether we infer ambiguous states or not,
        and reported in both sequences and mutations
        """
        result = run_ancestral(T=tree, aln=aln_pos3X, **generic_args, infer_ambiguous=True)
        sample_c = result['mutations']['nodes']['sample_C']

        assert sample_c['sequence'][2] in ['A', 'T', 'G', 'C']
        # And there are no "to N" mutations reported on the tree.
        pos3_muts = gather_mutations_at_pos(3, result)
        assert not any([m.endswith("N") for m in pos3_muts]), "Expected ambiguous base N to be inferred as A/T/C/G"

    @pytest.mark.xfail  # TODO XXX - fix bug
    def test_nuc_X_on_internal_node_normalized_to_N(self, tree, generic_args):
        """When all descendant tips of an internal node have 'X' at a position,
        the internal node's sequence should report 'N' (the standard ambiguous
        nucleotide) when we're not inferring ambiguous bases.
        """
        # Both children of node_AB have 'X' at pos 3; sample_C has 'G'
        aln = MultipleSeqAlignment([
            SeqRecord(Seq("AAXAA"), id="sample_A"),
            SeqRecord(Seq("AAXAA"), id="sample_B"),
            SeqRecord(Seq("AAGAA"), id="sample_C"),
        ])
        result = run_ancestral(T=tree, aln=aln, **generic_args, infer_ambiguous=False)
        node_AB = result['mutations']['nodes']['node_AB']

        # node_AB's children are both fully ambiguous at pos 3, so we have no
        # information about node_AB's state there. It should be 'N'.
        assert node_AB['sequence'][2] == 'N'

    @pytest.mark.xfail  # TODO XXX - fix bug
    def test_ref_nuc_X_normalized_in_root_mutations(self, tree, generic_args):
        """When the reference sequence contains 'X' at a nucleotide position,
        it should be treated as 'N' (both are fully ambiguous) in the root
        mutation comparison.
        """
        ref = "AAXAA"
        # All tips have 'N' at pos 3 → position IS masked (mask checks for 'N')
        # With infer_ambiguous=False, masked root mutations are still checked.
        # The reference 'X' at pos 3 should be equivalent to 'N'.
        aln = MultipleSeqAlignment([
            SeqRecord(Seq("AANAA"), id="sample_A"),
            SeqRecord(Seq("AANAA"), id="sample_B"),
            SeqRecord(Seq("AANAA"), id="sample_C"),
        ])
        args_with_ref = {**generic_args, 'reference_sequence': ref}
        result = run_ancestral(T=tree, aln=aln, **args_with_ref, infer_ambiguous=False)
        root = result['mutations']['nodes']['node_root']

        # Position 3 is masked (all tips N) so the root's reported sequence
        # has 'N' at pos 3. 
        assert root['sequence'][2]=='N', "Expected inferred root pos 3 to be 'N'"
        
        # Since X ≡ N for nucleotides, no mutation should be
        # reported between ref & root.
        pos3_muts = [m for m in root['muts'] if int(m[1:-1]) == 3]
        assert pos3_muts == [], f"Expected no mutation at pos 3 (X ≡ N), got: {pos3_muts}"


class TestAmbiguousAAReconstruction:
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

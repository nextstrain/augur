"""
Unit tests for nucleotide to acid translation
"""
from pathlib import Path
import sys

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation

# we assume (and assert) that this script is running from the tests/ directory
sys.path.append(str(Path(__file__).parent.parent.parent))

from augur import translate

class TestTranslate:
    def test_safe_translate(self):
        '''
        Test safe_translate from a given nucleotide sequence with no gaps to amino acid
        '''
        # initialize input and output tuples (tuple of params, expected outputs)
        params_and_outs = [(('ATG',), 'M'),
                           (('ATGGT-',), 'MX'),
                           (('ATG---',), 'M-'),
                           (('ATGTAG',), 'M*'),
                           (('',), ''),
                           (('ATGA-G',), 'MX')]

        # input each pair into the function and check
        for pair in params_and_outs:
            params, out = pair
            assert translate.safe_translate(*params) == out

    def test_safe_translate_errors(self):
        '''
        Test that safe_translate raises ValueError when sequence length is not divisible by 3
        '''
        invalid_sequences = [
            'A',
            'AT',
            'ATGT',
            'ATGTA',
            'ATGTAGA',
        ]

        for seq in invalid_sequences:
            with pytest.raises(ValueError, match="not divisible by 3"):
                translate.safe_translate(seq)

    def test_translate_feature(self):
        '''
        Test translate_feature from a dictionary of given nucleotides to dictionary of translated amino acids
        '''
        # Seq -> Amino https://en.wikipedia.org/wiki/DNA_codon_table
        seq1 = Seq("TTTCTTATGGTCGTA") 
        seq2 = Seq("TCTTCAACTGCTACA")
        seq3 = Seq("CATAATGAATATAAT")
        aln = {'seq1': seq1,
               'seq2': seq2,
               'seq3': seq3}
        feature = SeqFeature(FeatureLocation(0, 15), type="domain")

        # expected results
        expected_translations = {'seq1': 'FLMVV',
                                 'seq2': 'SSTAT',
                                 'seq3': 'HNEYN'}

        assert translate.translate_feature(aln, feature) == expected_translations

    # TODO: test_vcf_feature, assign_aa_vcf, assign_aa_fasta
    # Unclear how to emulate inputs (TreeTime dict, tree)


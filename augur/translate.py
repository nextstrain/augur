import os
import numpy as np
from Bio import SeqIO, SeqFeature, Seq, SeqRecord, Phylo
from .utils import read_node_data, load_features, write_json

def safe_translate(sequence, report_exceptions=False):
    """Returns an amino acid translation of the given nucleotide sequence accounting
    for gaps in the given sequence.

    Optionally, returns a tuple of the translated sequence and whether an
    exception was raised during initial translation.

    >>> safe_translate("ATG")
    'M'
    >>> safe_translate("ATGGT-")
    'MX'
    >>> safe_translate("ATG---")
    'M-'
    >>> safe_translate("ATGTAG")
    'M*'
    >>> safe_translate("")
    ''
    >>> safe_translate("ATGT")
    'M'
    >>> safe_translate("ATG", report_exceptions=True)
    ('M', False)
    >>> safe_translate("ATGA-G", report_exceptions=True)
    ('MX', True)
    """
    from Bio.Data.CodonTable import TranslationError
    from Bio.Seq import CodonTable
    translation_exception = False

    try:
        # Attempt translation by extracting the sequence according to the
        # BioPhython SeqFeature in frame gaps of three will translate as '-'
        translated_sequence = str(Seq.Seq(sequence).translate(gap='-'))
    except TranslationError:
        translation_exception = True

        # Any other codon like '-AA' or 'NNT' etc will fail. Translate codons
        # one by one.
        codon_table  = CodonTable.ambiguous_dna_by_name['Standard'].forward_table
        str_seq = str(sequence)
        codons = np.fromstring(str_seq[:len(str_seq) - len(str_seq) % 3], dtype='S3').astype('U')
        assert len(codons) > 0
        aas = []

        for c in codons:
            # Parse result of single codon translation, add amino acids as
            # appropriate.
            try:
                aa = codon_table.get(c)
                if aa is None:
                    if c == '---':
                        aas.append('-')
                    else:
                        aas.append('X')
                else:
                    aas.append(aa)
            except (TranslationError, ValueError):
                aas.append('X')

        translated_sequence = "".join(aas)

    if report_exceptions:
        return translated_sequence, translation_exception
    else:
        return translated_sequence


def translate_feature(aln, feature):
    '''
    translates a subsequence of input nucleotide sequences
    input:
        aln -- dictionary of sequences indexed by node name
        feature -- Biopython sequence feature
    returns:
        dictionary of translated sequences indexed by node name
    '''
    translations = {}
    for sname, seq in aln.items():
        aa_seq = safe_translate(str(feature.extract(seq)))
        translations[sname] = aa_seq

    return translations

def run(args):
    # read tree and data, if reading data fails, return with error code
    tree = Phylo.read(args.tree, 'newick')
    node_data = read_node_data(args.node_data)
    features = load_features(args.reference_sequence, args.genes)
    if features is None or node_data is None:
        print("ERROR: could not read node data (incl sequences) or features of reference sequence")
        return -1

    # extract sequences from node meta data
    sequences = {}
    for k,v in node_data['nodes'].items():
        if 'sequence' in v:
            sequences[k] = v['sequence']

    # translate every feature
    translations = {}
    for fname, feat in features.items():
        translations[fname] = translate_feature(sequences, feat)

    ## glob the annotations for later auspice export
    annotations = {}
    for fname, feat in features.items():
        annotations[fname] = {'start':int(feat.location.start),
                              'end':int(feat.location.end),
                              'strand': feat.location.strand}

    ## determine amino acid mutations for each node
    aa_muts = {}
    for n in tree.get_nonterminals():
        for c in n:
            aa_muts[c.name]={"aa_muts":{}}
        for fname, aln in translations.items():
            for c in n:
                if c.name in aln and n.name in aln:
                    tmp = [(a,pos,d) for pos, (a,d) in
                            enumerate(zip(aln[n.name], aln[c.name])) if a!=d]
                aa_muts[c.name]["aa_muts"][fname] = tmp

    write_json({'annotation':annotations, 'mutations':aa_muts}, args.output)

    ## write alignments to file is requested
    if args.alignment_output and '%GENE' in args.alignment_output:
        for fname, seqs in translations.items():
            SeqIO.write([SeqRecord.SeqRecord(seq=Seq.Seq(s), id=sname, name=sname, description='')
                         for sname, s in seqs.items()],
                         args.alignment_output.replace('%GENE', fname), 'fasta')


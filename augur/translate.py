"""
Translate gene regions from nucleotides to amino acids.
"""

import os, sys
import numpy as np
from Bio import SeqIO, SeqFeature, Seq, SeqRecord, Phylo
from .utils import read_node_data, load_features, write_json, write_VCF_translation, get_json_name
from treetime.vcf_utils import read_vcf

class MissingNodeError(Exception):
    pass

class MismatchNodeError(Exception):
    pass

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
    'MX'
    >>> safe_translate("ATG", report_exceptions=True)
    ('M', False)
    >>> safe_translate("ATGA-G", report_exceptions=True)
    ('MX', True)
    """
    from Bio.Data.CodonTable import TranslationError
    from Bio.Seq import CodonTable
    translation_exception = False

    #sequences not mod 3 give messy BiopythonWarning, so avoid by padding.
    if len(sequence)%3:
        sequence_padded = sequence + "N"*(3-len(sequence)%3)
    else:
        sequence_padded = sequence
    try:
        # Attempt translation by extracting the sequence according to the
        # BioPhython SeqFeature in frame gaps of three will translate as '-'
        translated_sequence = str(Seq.Seq(sequence_padded).translate(gap='-'))
    except TranslationError:
        translation_exception = True
        # Any other codon like '-AA' or 'NNT' etc will fail. Translate codons
        # one by one.
        codon_table  = CodonTable.ambiguous_dna_by_name['Standard'].forward_table
        str_seq = str(sequence_padded)
        codons = np.frombuffer(str_seq[:len(str_seq) - len(str_seq) % 3].encode(), dtype='S3').astype("U")
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
    Translates a subsequence of input nucleotide sequences.

    Parameters
    ----------
    aln : dict
        sequences indexed by node name

    feature : Bio.Seq.Seq
        BioPython sequence feature

    Returns
    -------
    dict :
        translated sequences indexed by node name

    '''
    translations = {}
    for sname, seq in aln.items():
        aa_seq = safe_translate(str(feature.extract(seq)))
        translations[sname] = aa_seq

    return translations


def translate_vcf_feature(sequences, ref, feature):
    '''Translates a subsequence of input nucleotide sequences.

    Parameters
    ----------
    sequences : dict
        TreeTime format dictionary from VCF-input of sequences indexed by node name

    ref :
        reference alignment the VCF was mapped to

    feature : Bio.Seq.Seq
        BioPython sequence feature

    Returns
    -------
    dict :
        translated reference gene, positions of AA differences, and AA
        differences indexed by node name

    '''
    def str_reverse_comp(str_seq):
        #gets reverse-compliment of a string and returns it as a string
        seq_str = Seq(str_seq)
        return str(seq_str.reverse_complement())

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    prot = {}
    prot['sequences'] = {}
    prot['positions'] = []

    refNuc = str(feature.extract( SeqRecord(seq=Seq(ref)) ).seq)
    # Need to get ref translation to store. check if multiple of 3 for sanity.
    # will be padded in safe_translate if not
    if len(refNuc)%3:
        print("Gene length of {} is not a multiple of 3. will pad with N".format(feature.qualifiers['Name'][0]), file=sys.stderr)

    ref_aa_seq = safe_translate(refNuc)
    prot['reference'] = ref_aa_seq

    start = int(feature.location.start) #1 less than GFF file; Bio 'pythons' the numbering
    end = int(feature.location.end)

    for seqk in sequences.keys():
        varSite = np.array(list(sequences[seqk].keys()))  #get positions where nuc diffs
        #reduce to only those within current gene
        geneVarSites = np.logical_and(varSite >= start, varSite < end)
        nucVarSites = varSite[geneVarSites] #translate this back to nuc position
        genNucSites = nucVarSites-start #get it in position within the gene to ensure in frame

        #Translate just the codon this nuc diff is in, and find out which AA loc
        #But need numbering to be w/in protin, not whole genome
        if feature.strand == -1:
            aaRepLocs = {(end-start-i-1)//3:safe_translate( str_reverse_comp( "".join([sequences[seqk][key+start]
                                    if key+start in sequences[seqk].keys() else ref[key+start]
                                for key in range(i-i%3,i+3-i%3)]) ))
                            for i in genNucSites}
        else:
            aaRepLocs = {i//3:safe_translate( "".join([sequences[seqk][key+start]
                                if key+start in sequences[seqk].keys() else ref[key+start]
                            for key in range(i-i%3,i+3-i%3)]) )
                        for i in genNucSites}

        aaRepLocsFinal = {}
        #remove if is a synonymous mutation!
        for key,val in aaRepLocs.items():
            if ref_aa_seq[key] != val:
                aaRepLocsFinal[key] = val
        aaRepLocs = aaRepLocsFinal

        #store the dict of differences
        prot['sequences'][seqk] = aaRepLocs
        #add to list of positions if needed
        for key in aaRepLocs.keys():
            if key not in prot['positions']:
                prot['positions'].append(key)

    prot['positions'].sort()

    #if no variable sites, exclude this gene
    if len(prot['positions']) == 0:
        return None
    else:
        return prot

def construct_mut(start, pos, end):
    return str(start) + str(pos) + str(end)

def assign_aa_vcf(tree, translations):
    aa_muts = {}

    #get mutations on the root
    root = tree.root
    aa_muts[root.name]={"aa_muts":{}}
    #If has no root node name, exit with error
    if root.name is None:
        print("\n*** Can't find node name for the tree root!")
        raise MissingNodeError()

    for fname, prot in translations.items():
        if root.name not in prot['sequences']:
            print("\n*** Can't find %s in the alignment provided!"%(root.name))
            raise MismatchNodeError()
        root_muts = prot['sequences'][root.name]
        tmp = []
        for pos in prot['positions']:
            if pos in root_muts:
                tmp.append(construct_mut(prot['reference'][pos], int(pos+1), root_muts[pos]))
        aa_muts[root.name]["aa_muts"][fname] = tmp

    for n in tree.get_nonterminals():
        for c in n:
            aa_muts[c.name]={"aa_muts":{}}
        for fname, prot in translations.items():
            if n.name not in prot['sequences']:
                print("\n*** Can't find %s in the alignment provided!"%(root.name))
                raise MismatchNodeError()
            n_muts = prot['sequences'][n.name]
            for c in n:
                tmp = []
                if c.name is None:
                    print("\n*** Internal node missing a node name!")
                    raise MissingNodeError()
                c_muts = prot['sequences'][c.name]
                for pos in prot['positions']:
                    #if pos in both, check if same
                    if pos in n_muts and pos in c_muts:
                        if n_muts[pos] != c_muts[pos]:
                            tmp.append(construct_mut(n_muts[pos], int(pos+1), c_muts[pos]))
                    elif pos in n_muts:
                        tmp.append(construct_mut(n_muts[pos], int(pos+1), prot['reference'][pos]))
                    elif pos in c_muts:
                        tmp.append(construct_mut(prot['reference'][pos], int(pos+1), c_muts[pos]))

                aa_muts[c.name]["aa_muts"][fname] = tmp

    return aa_muts

def assign_aa_fasta(tree, translations):
    aa_muts = {}

    #fasta input shouldn't have mutations on root, so give empty entry
    root = tree.get_nonterminals()[0]
    aa_muts[root.name]={"aa_muts":{}}

    for n in tree.get_nonterminals():
        if n.name is None:
            print("\n*** Tree is missing node names!")
            raise MissingNodeError()
        for c in n:
            aa_muts[c.name]={"aa_muts":{}}
        for fname, aln in translations.items():
            for c in n:
                if c.name in aln and n.name in aln:
                    tmp = [construct_mut(a, int(pos+1), d) for pos, (a,d) in
                            enumerate(zip(aln[n.name], aln[c.name])) if a!=d]
                    aa_muts[c.name]["aa_muts"][fname] = tmp
                elif c.name not in aln and n.name not in aln:
                    print("\n*** Can't find %s OR %s in the alignment provided!"%(c.name, n.name))
                    raise MismatchNodeError()
                else:
                    print("no sequence pair for nodes %s-%s"%(c.name, n.name))

        if n==tree.root:
            aa_muts[n.name]={"aa_muts":{}, "aa_sequences":{}}
            for fname, aln in translations.items():
                if n.name in aln:
                    aa_muts[n.name]["aa_sequences"][fname] = "".join(aln[n.name])

    return aa_muts

def get_genes_from_file(fname):
    genes = []
    if os.path.isfile(fname):
        with open(fname, encoding='utf-8') as ifile:
            for line in ifile:
                fields = line.strip().split('#')
                if fields[0].strip():
                    genes.append(fields[0].strip())
    else:
        print("File with genes not found. Looking for", fname)

    unique_genes = np.unique(np.array(genes))
    if len(unique_genes) != len(genes):
        print("You have duplicates in your genes file. They are being ignored.")
    print("Read in {} specified genes to translate.".format(len(unique_genes)))

    return unique_genes


def register_arguments(parser):
    parser.add_argument('--tree', help="prebuilt Newick -- no tree will be built if provided")
    parser.add_argument('--ancestral-sequences', type=str, help='JSON (fasta input) or VCF (VCF input) containing ancestral and tip sequences')
    parser.add_argument('--reference-sequence', required=True,
                        help='GenBank or GFF file containing the annotation')
    parser.add_argument('--genes', nargs='+', help="genes to translate (list or file containing list)")
    parser.add_argument('--output-node-data', type=str, help='name of JSON file to save aa-mutations to')
    parser.add_argument('--alignment-output', type=str, help="write out translated gene alignments. "
                                   "If a VCF-input, a .vcf or .vcf.gz will be output here (depending on file ending). If fasta-input, specify the file name "
                                   "like so: 'my_alignment_%%GENE.fasta', where '%%GENE' will be replaced by the name of the gene")
    parser.add_argument('--vcf-reference-output', type=str, help="fasta file where reference sequence translations for VCF input will be written")
    parser.add_argument('--vcf-reference', type=str, help='fasta file of the sequence the VCF was mapped to')


def run(args):
    ## read tree and data, if reading data fails, return with error code
    tree = Phylo.read(args.tree, 'newick')

    # If genes is a file, read in the genes to translate
    if args.genes and len(args.genes) == 1 and os.path.isfile(args.genes[0]):
        genes = get_genes_from_file(args.genes[0])
    else:
        genes = args.genes

    ## check file format and read in sequences
    is_vcf = False
    if any([args.ancestral_sequences.lower().endswith(x) for x in ['.vcf', '.vcf.gz']]):
        if not args.vcf_reference:
            print("ERROR: a reference Fasta is required with VCF-format input")
            return 1
        compress_seq = read_vcf(args.ancestral_sequences, args.vcf_reference)
        sequences = compress_seq['sequences']
        ref = compress_seq['reference']
        is_vcf = True
    else:
        node_data = read_node_data(args.ancestral_sequences, args.tree)
        if node_data is None:
            print("ERROR: could not read node data (incl sequences)")
            return 1
        # extract sequences from node meta data
        sequences = {}
        for k,v in node_data['nodes'].items():
            if 'sequence' in v:
                sequences[k] = v['sequence']

    ## load features; only requested features if genes given
    features = load_features(args.reference_sequence, genes)
    print("Read in {} features from reference sequence file".format(len(features)))
    if features is None:
        print("ERROR: could not read features of reference sequence file")
        return 1

    ### translate every feature - but not 'nuc'!
    translations = {}
    deleted = []
    for fname, feat in features.items():
        if is_vcf:
            trans = translate_vcf_feature(sequences, ref, feat)
            if trans:
                translations[fname] = trans
            else:
                deleted.append(fname)
        else:
            if feat.type != 'source':
                translations[fname] = translate_feature(sequences, feat)

    if len(deleted) != 0:
        print("{} genes had no mutations and so have been be excluded.".format(len(deleted)))

    ## glob the annotations for later auspice export
    #
    # Note that BioPython FeatureLocations use
    # "Pythonic" coordinates: [zero-origin, half-open)
    # Starting with augur v6 we use GFF coordinates: [one-origin, inclusive]
    annotations = {}
    for fname, feat in features.items():
        annotations[fname] = {'seqid':args.reference_sequence,
                              'type':feat.type,
                              'start':int(feat.location.start)+1,
                              'end':int(feat.location.end),
                              'strand': '+' if feat.location.strand else '-'}
    if is_vcf: #need to add our own nuc
        annotations['nuc'] = {'seqid':args.reference_sequence,
                              'type':feat.type,
                              'start': 1,
                              'end': len(ref),
                              'strand': '+'}

    ## determine amino acid mutations for each node
    try:
        if is_vcf:
            aa_muts = assign_aa_vcf(tree, translations)
        else:
            aa_muts = assign_aa_fasta(tree, translations)
    except MissingNodeError as err:
        print("\n*** ERROR: Some/all nodes have no node names!")
        print("*** Please check you are providing the tree output by 'augur refine'.")
        print("*** If you haven't run 'augur refine', please add node names to your tree by running:")
        print("*** augur refine --tree %s --output-tree <filename>.nwk"%(args.tree) )
        print("*** And use <filename>.nwk as the tree when running 'ancestral', 'translate', and 'traits'")
        return 1
    except MismatchNodeError as err:
        print("\n*** ERROR: Mismatch between node names in %s and in %s"%(args.tree, args.ancestral_sequences))
        print("*** Ensure you are using the same tree you used to run 'ancestral' as input here.")
        print("*** Or, re-run 'ancestral' using %s, then use the new %s as input here."%(args.tree, args.ancestral_sequences))
        return 1

    output_data = {'annotations':annotations, 'nodes':aa_muts}
    if is_vcf:
        output_data['reference'] = {}
        for fname in translations:
            output_data['reference'][fname] = translations[fname]['reference']
    else:
        output_data['reference'] = aa_muts[tree.root.name]['aa_sequences']

    out_name = get_json_name(args, '.'.join(args.tree.split('.')[:-1]) + '_aa-mutations.json')
    write_json(output_data, out_name)
    print("amino acid mutations written to", out_name, file=sys.stdout)

    ## write alignments to file is requested
    if args.alignment_output:
        if is_vcf:
            ## write VCF-style output if requested
            fileEndings = -1
            if args.alignment_output.lower().endswith('.gz'):
                fileEndings = -2
            vcf_out_ref = args.vcf_reference_output or '.'.join(args.alignment_output.split('.')[:fileEndings]) + '_reference.fasta'
            write_VCF_translation(translations, args.alignment_output, vcf_out_ref)
        else:
            ## write fasta-style output if requested
            if '%GENE' in args.alignment_output:
                for fname, seqs in translations.items():
                    SeqIO.write([SeqRecord.SeqRecord(seq=Seq.Seq(s), id=sname, name=sname, description='')
                                 for sname, s in seqs.items()],
                                 args.alignment_output.replace('%GENE', fname), 'fasta')
            else:
                print("ERROR: alignment output file does not contain '%GENE', so will not be written.")

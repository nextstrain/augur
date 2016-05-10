'''
parse, subsample, and align a sequence data set
'''
from __future__ import division, print_function
import os, re, time, csv, sys
from io_util import myopen, make_dir, remove_dir, tree_to_json, write_json
from collections import defaultdict
from Bio import SeqIO
import numpy as np
from seq_util import pad_nucleotide_sequences, nuc_alpha, aa_alpha
from datetime import datetime

TINY = 1e-10

def fix_names(n):
    return n.replace(" ","_").replace("(",'_').replace(")",'_').replace("'",'_').replace(":",'_')

def calc_af(aln, alpha):
    aln_array = np.array(aln)
    af = np.zeros((len(alpha), aln_array.shape[1]))
    for ai, state in enumerate(alpha):
        af[ai] += (aln_array==state).mean(axis=0)
    af[-1] = 1.0 - af[:-1].sum(axis=0)
    return af

def num_date(date):
    days_in_year = date.toordinal()- datetime(year=date.year, month=1, day=1).date().toordinal()
    return date.year + days_in_year/365.25

class sequence_set(object):
    """sequence_set subsamples a set of sequences, aligns them and exports variability statistics"""
    def __init__(self, fname=None, reference= None, **kwarks):
        super(sequence_set, self).__init__()
        self.kwarks = kwarks
        self.nthreads = 2
        if fname is not None and os.path.isfile(fname):
            with myopen(fname) as seq_file:
                self.raw_seqs = {fix_names(x.description):x for x in SeqIO.parse(seq_file, 'fasta')}
                for x in self.raw_seqs.values():
                    x.id = fix_names(x.id)
                    x.name = fix_names(x.id)
                    x.description = fix_names(x.description)
        if 'run_dir' not in kwarks:
            import random
            self.run_dir = '_'.join(['temp', time.strftime('%Y%m%d-%H%M%S',time.gmtime()), str(random.randint(0,1000000))])
        else:
            self.run_dir = kwarks['run_dir']

        if reference is not None:
            if type(reference) is str and fix_names(reference) in self.raw_seqs:
                self.reference = self.raw_seqs[fix_names(reference)]
            else:
                self.reference = reference
        else: self.reference=None

    def parse(self, fields, sep='|', strip='_'):
        '''
        split the sequence description and add annotations to sequences
        '''
        for seq in self.raw_seqs.values():
            if not hasattr(seq, "attributes"): seq.attributes = {}
            words = map(lambda x:x.strip(strip),seq.description.replace(">","").split(sep))
            for ii, val in enumerate(words):
                if ii in fields:
                    if val not in ["", "-"]:
                        seq.attributes[fields[ii]] = val
                    else:
                        seq.attributes[fields[ii]] = ""

    def from_vdb(self, virus):
        from nextstraindb.vdb.src.vdb_download import vdb_download
        myvdb = vdb_download(virus=virus)
        myvdb.vdb_download(output=False)
        for v in myvdb.viruses:
            attr = {v[k] for k in v if k!='sequence'}
            seq = SeqRecord(Seq(v['sequence']), id=attr['strain'], name=attr['strain'], description=attr['strain'])
            seq.attributes = attr
            self.raw_seqs[attr['strain']] = seq

    def ungap(self):
        '''
        remove previously existing gaps and make sure all sequences are upper case
        '''
        for seq in self.raw_seqs.values():
            seq.seq = seq.seq.ungap('-').upper()

    def parse_date(self, fmts, prune=True):
        if not hasattr(self.raw_seqs.values()[0], "attributes"):
            print("parse meta info first")
            return
        from datetime import datetime
        for seq in self.raw_seqs.values():
            if 'date' in seq.attributes and seq.attributes['date']!='':
                for fmt in fmts:
                    try:
                        if callable(fmt):
                            tmp = fmt(seq.attributes['date'])
                        else:
                            tmp = datetime.strptime(seq.attributes['date'], fmt).date()
                        seq.attributes['raw_date'] = seq.attributes['date']
                        seq.attributes['num_date'] = num_date(tmp)
                        seq.attributes['date']=tmp
                        break
                    except:
                        continue

        if prune:
            self.raw_seqs = {k:v for k,v in self.raw_seqs.iteritems()
                            if 'date' in v.attributes and type(v.attributes['date'])!=str}

    def filter(self, func):
        self.raw_seqs = {key:seq for key, seq in self.raw_seqs.iteritems() if func(seq)}

    def clock_filter(self, root_seq=None, n_iqd=3, max_gaps = 1.0, plot=False):
        '''
        remove sequences form the set that are that evolve much faster or slower
        compared the majority. Regions with predominantly gaps can be removed since
        this can skew the evolutionary rates.
        '''
        from Bio.Align import MultipleSeqAlignment
        if root_seq is None: # use consensus
            af = calc_af(self.aln, nuc_alpha)
            root_seq = np.fromstring(nuc_alpha, 'S1')[af.argmax(axis=0)]
        if type(root_seq)==str and root_seq in self.sequence_lookup:
            root_seq = np.array(self.sequence_lookup[root_seq])
        if max_gaps<1.0:
            af=calc_af(self.aln, nuc_alpha)
            good_pos = af[nuc_alpha.index('-')]<max_gaps
        else:
            good_pos = np.ones(self.aln.get_alignment_length(), dtype=bool)
        date_vs_distance = {}
        for seq in self.aln:
            date_vs_distance[seq.id] = (seq.attributes['num_date'],
                np.mean((np.array(seq)!=root_seq)[(np.array(seq)!='-')&(root_seq!='-')&good_pos]))
        date_vs_distance_array=np.array(date_vs_distance.values())
        from scipy.stats import linregress, scoreatpercentile
        slope, intercept, rval, pval, stderr = linregress(date_vs_distance_array[:,0], date_vs_distance_array[:,1])
        print("distance vs time regression:",slope)
        residuals = (intercept + slope*date_vs_distance_array[:,0]) - date_vs_distance_array[:,1]
        IQD = scoreatpercentile(residuals, 75) - scoreatpercentile(residuals,25)
        if plot:
            import matplotlib.pyplot as plt
            plt.ion()
            plt.scatter(date_vs_distance_array[:,0], date_vs_distance_array[:,1], c='g')
            bad_points = abs(intercept+slope*date_vs_distance_array[:,0] - date_vs_distance_array[:,1])>n_iqd*IQD
            plt.scatter(date_vs_distance_array[bad_points,0], date_vs_distance_array[bad_points,1], c='r')


        print("before clock filter:",len(self.aln))
        self.aln = MultipleSeqAlignment([seq for seq in self.aln
                    if abs(intercept+slope*date_vs_distance[seq.id][0] - date_vs_distance[seq.id][1])<n_iqd*IQD])
        print("after clock filter:",len(self.aln))

    def subsample(self, category=None, priority=None, threshold=None, repeated=False):
        '''
        produce a useful set of sequences from the raw input.
        arguments:
        category  -- callable that assigns each sequence to a category for subsampling
        priority  -- callable that assigns each sequence a priority to be included in
                     the final sample. this is applied independently in each category
        threshold -- callable that determines the number of sequences from each category
                     that is included in the final set. takes arguments, cat and seq
                     alternatively can be an int
        '''
        if category is None:
            category = lambda x:(x.attributes['date'].year, x.attributes['date'].month)
        if priority is None:
            priority = lambda x:np.random.random()
        if threshold is None:
            threshold = lambda x:5
        elif type(threshold) is int:
            print("using threshold:",threshold)
            tmp = threshold
            threshold = lambda x:tmp

        self.sequence_categories = defaultdict(list)
        if repeated:
            seqs_to_subsample = self.seqs.values()
        else:
            seqs_to_subsample = self.raw_seqs.values()

        for seq in seqs_to_subsample:
            seq._priority = priority(seq)
            self.sequence_categories[category(seq)].append(seq)

        self.seqs = {}
        for cat, seqs in self.sequence_categories.iteritems():
            under_sampling = min(1.00, 1.0*len(seqs)/threshold)
            for s in seqs: s.under_sampling=under_sampling
            seqs.sort(key=lambda x:x._priority, reverse=True)
            self.seqs.update({seq.id:seq for seq in seqs[:threshold( (cat, seqs) )]})

        if self.reference.id not in self.seqs:
            self.seqs[self.reference.id] = self.reference

    def align(self):
        from Bio import AlignIO
        make_dir(self.run_dir)
        os.chdir(self.run_dir)

        SeqIO.write(self.seqs.values(), "temp_in.fasta", "fasta")
        os.system("mafft --anysymbol --thread " + str(self.nthreads) + " temp_in.fasta > temp_out.fasta")

        self.aln = AlignIO.read('temp_out.fasta', 'fasta')
        self.sequence_lookup = {seq.id:seq for seq in self.aln}
        self.reference_aligned = self.sequence_lookup[self.reference.id]
        # add attributes to alignment
        for seqid, seq in self.seqs.iteritems():
            self.sequence_lookup[seqid].attributes = seq.attributes
        os.chdir('..')
        remove_dir(self.run_dir)

    def codon_align(self, alignment_tool="mafft", prune=True):
        ''' takes a nucleotide alignment, translates it, aligns the amino acids, pads the gaps
        note that this suppresses any compensated frameshift mutations

        Parameters:
        - alignment_tool: ['mafft', 'muscle'] the commandline tool to use
        '''
        from Bio import AlignIO,SeqIO
        from Bio.SeqRecord import SeqRecord
        make_dir(self.run_dir)
        os.chdir(self.run_dir)

        # translage
        aa_seqs = {}
        for seq in self.seqs.values():
            tempseq = seq.seq.translate()
            # use only sequences that translate with out trouble
            if '*' not in str(tempseq)[:-1] or prune==False:
                aa_seqs[seq.id]=SeqRecord(tempseq,id=seq.id)
            else:
                print(seq.id,"has premature stops, discarding")

        tmpfname = 'temp_in.fasta'
        SeqIO.write(aa_seqs.values(), tmpfname,'fasta')

        if alignment_tool=='muscle':
            from Bio.Align.Applications import MuscleCommandline
            cline = MuscleCommandline(input=tmpfname, out=tmpfname[:-5]+'aligned.fasta')
            cline()
            aln_aa = AlignIO.read(tmpfname[:-5]+'aligned.fasta', "fasta")
        elif alignment_tool=='mafft':
            from Bio.Align.Applications import MafftCommandline
            from StringIO import StringIO
            mafft_cline = MafftCommandline(input=tmpfname)
            stdout, stderr = mafft_cline()
            aln_aa = AlignIO.read(StringIO(stdout), "fasta")
        else:
            print('Alignment tool not supported:',alignment_tool)
            return

        #generate nucleotide alignment
        self.aln = pad_nucleotide_sequences(aln_aa, self.seqs)
        self.sequence_lookup = {seq.id:seq for seq in self.aln}
        self.reference_aligned = self.sequence_lookup[self.reference.id]
        # add attributes to alignment
        for seq in self.seqs.values():
            if seq.id in self.sequence_lookup:
                self.sequence_lookup[seq.id].attributes = seq.attributes
        os.chdir('..')
        remove_dir(self.run_dir)


    def strip_non_reference(self):
        ungapped = np.array(self.reference_aligned)!='-'
        from Bio.Seq import Seq
        for seq in self.aln:
            seq.seq = Seq("".join(np.array(seq)[ungapped]))

    def diversity_statistics(self):
        ''' calculate alignment entropy of nucleotide and optionally protein alignments '''
        if not hasattr(self, "aln"):
            print("calculate alignment first")
            return
        aln_array = np.array(self.aln)
        self.af = {'nuc': calc_af(self.aln, nuc_alpha)}
        tmp_af = self.af['nuc'][:-2]/self.af['nuc'][:-2].sum(axis=0)
        self.entropy ={'nuc': -(tmp_af*np.log(tmp_af+TINY)).sum(axis=0)}

        if hasattr(self, "translations"):
            for prot, aln in self.translations.iteritems():
                self.af[prot] = calc_af(aln, aa_alpha)
                tmp_af = self.af[prot][:-2]/self.af[prot][:-2].sum(axis=0)
                self.entropy[prot] = -(tmp_af*np.log(tmp_af+TINY)).sum(axis=0)

    def translate(self, proteins=None):
        from Bio.SeqFeature import FeatureLocation
        from Bio.Seq import Seq
        from Bio.Align import MultipleSeqAlignment
        if not hasattr(self, "proteins"): # generate dictionaries to hold annotation and translation
            self.translations={}
            self.proteins={}

        if proteins is None: # add a default translation of the entire sequence unless otherwise specified
            self.proteins.update({'cds':FeatureLocation(start=0, end=self.aln.get_alignment_length(), strand=1)})
        else:
            self.proteins.update(proteins)

        for prot in self.proteins:
            aa_seqs = []
            for seq in self.aln:
                try:
                    # soon not needed as future biopython version will translate --- into -
                    tmpseq = self.proteins[prot].extract(seq)
                    tmpseq.attributes = seq.attributes
                    tmpseq.seq = Seq(str(Seq(str(tmpseq.seq).replace('---', 'NNN')).translate()).replace('X','-'))
                except:
                    tmpseq.seq = Seq(str(Seq("".join([x if x in 'ACGT' else 'N' for x in str(tmpseq.seq)])).translate()).replace('X','-'))
                    print("Trouble translating",seq.id)
                    #import ipdb; ipdb.set_trace()
                aa_seqs.append(tmpseq)
            self.translations[prot] = MultipleSeqAlignment(aa_seqs)

    def export_diversity(self, fname = 'entropy.json'):
        if not hasattr(self, "entropy"):
            self.diversity_statistics()
        entropy_json = {}
        for feat in self.entropy:
            S = [max(0,round(x,4)) for x in self.entropy[feat]]
            n = len(S)
            if feat=='nuc':
                entropy_json[feat] = {'pos':range(0,n), 'codon':[x//3 for x in range(0,n)], 'val':S}
            else:
                entropy_json[feat] = {'pos':[x for x in self.proteins[feat]][::3],
                                      'codon':[(x-self.proteins[feat].start)//3 for x in self.proteins[feat]][::3], 'val':S}
        write_json(entropy_json, fname, indent=None)

if __name__=="__main__":
    myseqs = sequence_set('data/gag.fasta.gz', reference='B|FR|1985|NL4_3_LAI_NY5_pNL43_NL43|244167|NL43|325|U26942')
    myseqs.ungap()
    myseqs.parse({0:"subtype", 1:"country", 2:"date", 4:"name", 5:"id", 6:"patient", 7:"accession"})
    myseqs.parse_date(["%Y-%m-%d", "%Y"])
    myseqs.filter(lambda x:x.attributes['subtype']=='B')
    myseqs.subsample(category = lambda x:x.attributes['date'].year)
    myseqs.codon_align(prune=True)
    #myseqs.align()
    myseqs.strip_non_reference()
    myseqs.clock_filter(n_iqd=3, plot=True)
    myseqs.diversity_statistics()
    myseqs.translate()



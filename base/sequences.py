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

def ambiguous_date_to_date_range(mydate, fmt):
    sep = fmt.split('%')[1][-1]
    min_date, max_date = {}, {}
    for val, field  in zip(mydate.split(sep), fmt.split(sep+'%')):
        f = 'year' if 'y' in field.lower() else ('day' if 'd' in field.lower() else 'month')
        if 'XX' in val:
            if f=='year':
                return None, None
            elif f=='month':
                min_date[f]=1
                max_date[f]=12
            elif f=='day':
                min_date[f]=1
                max_date[f]=31
        else:
            min_date[f]=int(val)
            max_date[f]=int(val)
    max_date['day'] = min(max_date['day'], 31 if max_date['month'] in [1,3,5,7,8,10,12]
                                           else 28 if max_date['month']==2 else 30)
    return (datetime(year=min_date['year'], month=min_date['month'], day=min_date['day']).date(),
            datetime(year=max_date['year'], month=max_date['month'], day=max_date['day']).date())

class sequence_set(object):
    """sequence_set subsamples a set of sequences, aligns them and exports variability statistics"""
    def __init__(self, fname=None, reference_seq= None, **kwarks):
        super(sequence_set, self).__init__()
        self.kwarks = kwarks
        self.nthreads = 2
        if fname is not None and os.path.isfile(fname):
            with myopen(fname) as seq_file:
                self.all_seqs = {x.name:x for x in SeqIO.parse(seq_file, 'fasta')}
        elif 'virus' in kwarks:
            self.from_vdb(kwarks['virus'])
        else:
            print('no input sequences found -- empty sequence set')
            return

        if 'run_dir' not in kwarks:
            import random
            self.run_dir = '_'.join(['temp', time.strftime('%Y%m%d-%H%M%S',time.gmtime()),
                                     str(random.randint(0,1000000))])
        else:
            self.run_dir = kwarks['run_dir']

        if reference_seq is not None:
            if type(reference_seq) is str and reference_seq in self.all_seqs:
                self.reference_seq = self.all_seqs[reference_seq]
            else:
                self.reference_seq = reference_seq
        else:
            self.reference_seq=None


    def parse(self, fields, sep='|', strip='_'):
        '''
        split the sequence description and add annotations to sequences
        '''
        for seq in self.all_seqs.values():
            if not hasattr(seq, "attributes"): seq.attributes = {}
            words = map(lambda x:x.strip(strip),seq.description.replace(">","").split(sep))
            for ii, val in enumerate(words):
                if ii in fields:
                    if val not in ["", "-"]:
                        seq.attributes[fields[ii]] = val
                    else:
                        seq.attributes[fields[ii]] = ""
        if 'strain' in fields.values():
            self.all_seqs = {seq.attributes['strain']:seq for seq in self.all_seqs.values()}
            for seq in self.all_seqs.values():
                seq.id = seq.attributes['strain']
                seq.name = seq.attributes['strain']


    def ungap(self):
        '''
        remove previously existing gaps and make sure all sequences are upper case
        '''
        for seq in self.all_seqs.values():
            seq.seq = seq.seq.ungap('-').upper()

    def parse_date(self, fmts, prune=True):
        if not hasattr(self.all_seqs.values()[0], "attributes"):
            print("parse meta info first")
            return
        from datetime import datetime
        for seq in self.all_seqs.values():
            if 'date' in seq.attributes and seq.attributes['date']!='':
                for fmt in fmts:
                    try:
                        if 'XX' in seq.attributes['date']:
                            min_date, max_date = ambiguous_date_to_date_range(seq.attributes['date'], fmt)
                            seq.attributes['raw_date'] = seq.attributes['date']
                            seq.attributes['num_date'] = np.array((num_date(min_date), num_date(max_date)))
                            seq.attributes['date'] = min_date
                        else:
                            if callable(fmt):
                                tmp = fmt(seq.attributes['date'])
                            else:
                                try:
                                    tmp = datetime.strptime(seq.attributes['date'], fmt).date()
                                except:
                                    tmp = seq.attributes['date']
                            seq.attributes['raw_date'] = seq.attributes['date']
                            seq.attributes['num_date'] = num_date(tmp)
                            seq.attributes['date']=tmp
                            break
                    except:
                        continue

        if prune:
            self.filter(func = lambda x:'date' in x.attributes and type(x.attributes['date'])!=str)


    def filter(self, func, leave_ref=False):
        if leave_ref:
            self.all_seqs = {key:seq for key, seq in self.all_seqs.iteritems() if func(seq) or key==self.reference_seq.name}
        else:
            self.all_seqs = {key:seq for key, seq in self.all_seqs.iteritems() if func(seq)}
        print("Filtered to %d sequences"%len(self.all_seqs))

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
        self.reference_aln = None
        for seq in self.aln:
            date_vs_distance[seq.id] = (seq.attributes['num_date'],
                np.mean((np.array(seq)!=root_seq)[(np.array(seq)!='-')&(root_seq!='-')&good_pos]))
            if seq.id==self.reference.id:
                self.reference_aln = seq
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
        tmp = {seq.id:seq for seq in self.aln
                if abs(intercept+slope*date_vs_distance[seq.id][0] - date_vs_distance[seq.id][1])<n_iqd*IQD}
        if self.reference.id not in tmp and self.reference_aln is not None:
            print('adding reference again after clock filter')
            tmp[self.reference.id] = self.reference_aln
        self.aln = MultipleSeqAlignment(tmp.values())
        print("after clock filter:",len(self.aln))


    def subsample(self, category=None, priority=None, threshold=None, repeated=False, forced_strains=[]):
        '''
        produce a useful set of sequences from the raw input.
        arguments:
        category  -- callable that assigns each sequence to a category for subsampling
        priority  -- callable that assigns each sequence a priority to be included in
                     the final sample. this is applied independently in each category
        threshold -- callable that determines the number of sequences from each category
                     that is included in the final set. takes arguments, cat and seq
                     alternatively can be an int
        forced_strains -- list of of strain names that should always be included (set to high priorty)
        '''
        # define filter criteria if not specified
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

        # if we do repeated subsampling, subsamples seqs, otherwise all_seqs
        self.sequence_categories = defaultdict(list)
        if repeated:
            seqs_to_subsample = self.seqs.values()
        else:
            seqs_to_subsample = self.all_seqs.values()

        # sort sequences into categories and assign priority score
        for seq in seqs_to_subsample:
            seq._priority = priority(seq)
            if seq.id in forced_strains:
                seq._priority = 1.0
            self.sequence_categories[category(seq)].append(seq)

        # sample and record the degree to which a category is under_sampled
        self.seqs = {}
        for cat, seqs in self.sequence_categories.iteritems():
            under_sampling = min(1.00, 1.0*len(seqs)/threshold(cat))
            for s in seqs: s.under_sampling=under_sampling
            seqs.sort(key=lambda x:x._priority, reverse=True)
            self.seqs.update({seq.id:seq for seq in seqs[:threshold( (cat, seqs) )]})

        print("Subsampled to %d sequences"%len(self.all_seqs))
        print("Subsampled to %d sequences"%len(self.seqs))

    def align(self, debug=False):
        '''
        align sequences using mafft
        '''
        from Bio import AlignIO
        from Bio.Align import MultipleSeqAlignment
        make_dir(self.run_dir)
        os.chdir(self.run_dir)
        ref_in_set = self.reference_seq.name in self.seqs
        if ref_in_set:
            out_seqs = self.seqs.values()
        else:
            out_seqs = self.seqs.values() + [self.reference_seq]
        print("align: reference in set",ref_in_set)
        SeqIO.write(out_seqs, "temp_in.fasta", "fasta")
        os.system("mafft --anysymbol --thread " + str(self.nthreads) + " temp_in.fasta 1> temp_out.fasta 2>mafft_stderr")

        tmp_aln = AlignIO.read('temp_out.fasta', 'fasta')
        self.sequence_lookup = {seq.id:seq for seq in tmp_aln}
        # add attributes to alignment
        for seqid, seq in self.seqs.iteritems():
            self.sequence_lookup[seqid].attributes = seq.attributes
        self.aln = MultipleSeqAlignment([s for s in tmp_aln
                            if s.name!=self.reference_seq.name or ref_in_set])
        os.chdir('..')
        if not debug:
            remove_dir(self.run_dir)


    def codon_align(self, alignment_tool="mafft", prune=True, verbose=0, debug=False):
        ''' takes a nucleotide alignment, translates it, aligns the amino acids, pads the gaps
        note that this suppresses any compensated frameshift mutations

        Parameters:
        - alignment_tool: ['mafft', 'muscle'] the commandline tool to use
        '''
        from Bio import AlignIO,SeqIO
        from Bio.SeqRecord import SeqRecord
        make_dir(self.run_dir)
        os.chdir(self.run_dir)

        # translate
        aa_seqs = {}
        bad_seq = 0
        for seq in self.seqs.values():
            tempseq = seq.seq.translate()
            # use only sequences that translate with out trouble
            if '*' not in str(tempseq)[:-1] or prune==False:
                aa_seqs[seq.id]=SeqRecord(tempseq,id=seq.id)
                aa_seqs[seq.id].attributes = seq.attributes
            else:
                if verbose: print(seq.id,"has premature stops, discarding")
            bad_seq+='*' in str(tempseq)[:-1]

        print('Number of sequences with stops:',bad_seq,'out of total',len(self.seqs))
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
        # add attributes to alignment
        for seq in self.seqs.values():
            if seq.id in self.sequence_lookup:
                self.sequence_lookup[seq.id].attributes = seq.attributes
        os.chdir('..')
        if not debug:
            remove_dir(self.run_dir)


    def strip_non_reference(self):
        ungapped = np.array(self.sequence_lookup[self.reference_seq.name])!='-'
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
        '''
        make alignment of translations
        '''
        from Bio.SeqFeature import FeatureLocation
        from Bio.Seq import Seq
        from Bio.Align import MultipleSeqAlignment
        if not hasattr(self, "proteins"): # generate dictionaries to hold annotation and translation
            self.translations={}
            self.proteins={}

        # add a default translation of the entire sequence unless otherwise specified
        if proteins is None and len(self.proteins)==0:
            self.proteins.update({'cds':FeatureLocation(start=0,
                end=self.aln.get_alignment_length(), strand=1)})
        else:
            self.proteins.update(proteins)

        for prot in self.proteins:
            aa_seqs = []
            for seq in self.aln:
                try:
                    # soon not needed as future biopython version will translate --- into -
                    tmpseq = self.proteins[prot].extract(seq)
                    tmpseq.attributes = seq.attributes
                    tmpseq.seq = Seq(str(Seq(str(tmpseq.seq).replace('---', 'NNN'))
                                         .translate()).replace('X','-'))
                except:
                    tmpseq.seq = Seq(str(Seq("".join([x if x in 'ACGT' else 'N'
                        for x in str(tmpseq.seq)])).translate()).replace('X','-'))
                    print("Trouble translating",seq.id)
                aa_seqs.append(tmpseq)
            self.translations[prot] = MultipleSeqAlignment(aa_seqs)


    def export_diversity(self, fname = 'entropy.json', indent=None):
        '''
        write the alignment entropy of each alignment (nucleotide and translations) to file
        '''
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
        write_json(entropy_json, fname, indent=indent)


if __name__=="__main__":
    pass

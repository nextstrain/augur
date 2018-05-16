from Bio import SeqIO
import pandas as pd
from collections import defaultdict
import random,os
from .utils import read_metadata, get_numerical_dates

comment_char = '#'

def run(args):
    '''
    filter and subsample a set of sequences into an analysis set
    '''
    seqs = {seq.id:seq for seq in SeqIO.parse(args.sequences, 'fasta')}
    seq_names = list(seqs.keys())
    meta_dict, meta_columns = read_metadata(args.metadata)

    if args.exclude and os.path.isfile(args.exclude):
        with open(args.exclude, 'r') as ifile:
            to_exclude = set([line.strip() for line in ifile if line[0]!=comment_char])
        seq_names = [s for s in seq_names if s not in to_exclude]

    if (args.min_date or args.max_date) and 'date' in meta_columns:
        dates = get_numerical_dates(metadata)
        if args.min_date:
            seq_names = [s for s in seq_names if np.min(dates[s])>args.min_date]
        if args.max_date:
            seq_names = [s for s in seq_names if np.min(dates[s])<args.max_date]

    if args.cat and args.viruses_per_cat:
        vpc = args.viruses_per_cat
        seq_names_by_cat = defaultdict(list)

        for seq_name in seq_names:
            cat = []
            if seq_name not in meta_dict:
                print("WARNING: no metadata for %s, skipping"%seq_name)
                continue
            else:
                m = meta_dict[seq_name]
            for c in args.cat:
                if c in m:
                    cat.append(m.loc[c])
                elif c in ['month', 'year'] and 'date' in m:
                    try:
                        year = int(m["date"].split('-')[0])
                    except:
                        print("WARNING: no valid year, skipping",seq_name, m["date"])
                        continue
                    if c=='month':
                        try:
                            month = int(m["date"].split('-')[1])
                        except:
                            month = random.randint(1,12)
                        cat.append((year, month))
                    else:
                        cat.append(year)
            seq_names_by_cat[tuple(cat)].append(seq_name)

        seq_subsample = []
        for cat, s in seq_names_by_cat.items():
            tmp_seqs = [seqs[seq_name] for seq_name in s]
            seq_subsample.extend(tmp_seqs if len(s)<=vpc
                                 else random.sample(tmp_seqs,vpc))
    else:
        seq_subsample = list(seqs.values())

    SeqIO.write(seq_subsample, args.output, 'fasta')



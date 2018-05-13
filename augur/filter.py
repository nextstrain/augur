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
    seqs = list(SeqIO.parse(args.sequences, 'fasta'))
    meta_dict, meta_columns = read_metadata(args.metadata)

    if args.exclude and os.path.isfile(args.exclude):
        with open(args.exclude, 'r') as ifile:
            to_exclude = set([line.strip() for line in ifile if line[0]!=comment_char])
        seqs = [s for s in seqs if s.id not in to_exclude]

    if (args.min_date or args.max_date) and 'date' in meta_columns:
        dates = get_numerical_dates(metadata)
        if args.min_date:
            seqs = [s for s in seq if np.min(dates[s.id])>args.min_date]
        if args.max_date:
            seqs = [s for s in seq if np.min(dates[s.id])<args.max_date]

    if args.cat and args.viruses_per_cat:
        vpc = args.viruses_per_cat
        seqs_by_cat = defaultdict(list)

        for seq in seqs:
            cat = []
            if seq.id not in meta_dict:
                print("WARNING: no metadata for %s, skipping"%seq.id)
                continue
            else:
                m = meta_dict[seq.id]
            for c in args.cat:
                if c in m:
                    cat.append(m.loc[c])
                elif c in ['month', 'year'] and 'date' in m:
                    try:
                        year = int(m["date"].split('-')[0])
                    except:
                        print("WARNING: no valid year, skipping",seq.id, m["date"])
                        continue
                    if c=='month':
                        try:
                            month = int(m["date"].split('-')[1])
                        except:
                            month = random.randint(1,12)
                        cat.append((year, month))
                    else:
                        cat.append(year)
            seqs_by_cat[tuple(cat)].append(seq)

        seqs = []
        for cat, s in seqs_by_cat.items():
            seqs.extend(s if len(s)<=vpc else random.sample(s,vpc))

    SeqIO.write(seqs, args.output, 'fasta')



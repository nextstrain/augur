from __future__ import division, print_function
import matplotlib as mpl
mpl.use('Agg')
from collections import defaultdict
import sys
sys.path.insert(0,'.')  # need to import from base
print(sys.argv, sys.path)
from base.process import process
import numpy as np
from datetime import datetime, timedelta
from base.io_util import myopen

regions = ['africa', 'south_asia', 'europe', 'china', 'north_america',
           'china', 'south_america', 'japan_korea', 'oceania', 'southeast_asia']

region_cmap = [
    ["africa",          "#5097BA"],
    ["south_america",   "#60AA9E"],
    ["west_asia",       "#75B681"],
    ["oceania",         "#8EBC66"],
    ["europe",          "#AABD52"],
    ["japan_korea",     "#C4B945"],
    ["north_america",   "#D9AD3D"],
    ["southeast_asia",  "#E59637"],
    ["south_asia",      "#E67030"],
    ["china",           "#DF4327"]
];

region_groups = {'NA':'north_america',
                 'AS':['china', 'japan_korea', 'south_asia', 'southeast_asia'],
                 'OC':'oceania', 'EU':'europe'}

attribute_nesting = {'geographic location':['region', 'country', 'city'],}

geo_attributes = ['region']


color_options = {
    "country":{"key":"country", "legendTitle":"Country", "menuItem":"country", "type":"discrete"},
    "region":{"key":"region", "legendTitle":"Region", "menuItem":"region", "type":"discrete", "color_map":region_cmap},
    "num_date":{"key":"num_date", "legendTitle":"Sampling date", "menuItem":"date", "type":"continuous"},
    "ep":{"key":"ep", "legendTitle":"Epitope Mutations", "menuItem":"epitope mutations", "type":"continuous"},
    "ne":{"key":"ne", "legendTitle":"Non-epitope Mutations", "menuItem":"nonepitope mutations", "type":"continuous"},
    "rb":{"key":"rb", "legendTitle":"Receptor Binding Mutations", "menuItem":"RBS mutations", "type":"continuous"},
    "gt":{"key":"genotype", "legendTitle":"Genotype", "menuItem":"genotype", "type":"discrete"},
    "cHI":{"key":"cHI", "legendTitle":"Antigenic advance", "menuItem":"Antigenic", "type":"continuous"}
}
panels = ['tree', 'entropy', 'frequencies']

outliers = {
'h3n2':["A/Sari/388/2006", "A/SaoPaulo/36178/2015", "A/Pennsylvania/40/2010", "A/Pennsylvania/14/2010",
        "A/Pennsylvania/09/2011", "A/OSAKA/31/2005", "A/Ohio/34/2012", "A/Kenya/170/2011", "A/Kenya/168/2011",
        "A/Indiana/21/2013", "A/Indiana/13/2012", "A/Indiana/6/2013", "A/Indiana/17/2013", "A/Indiana/11/2013", "A/Indiana/08/2012", "A/Indiana/06/2013",
        "A/India/6352/2012", "A/HuNan/01/2014", "A/Helsinki/942/2013", "A/Guam/AF2771/2011", "A/Chile/8266/2003",
        "A/Busan/15453/2009", "A/Nepal/142/2011", "A/Kenya/155/2011", "A/Guam/AF2771/2011", "A/Michigan/82/2016",
        "A/Ohio/27/2016", "A/Ohio/28/2016", "A/Michigan/83/2016", "A/Michigan/84/2016", "A/Jiangsu-Tianning/1707/2013",
        "A/HuNan/1/2014", "A/Iran/227/2014", "A/Iran/234/2014", "A/Iran/140/2014", "A/Jiangsu-Chongchuan/1830/2014",
        "A/Chile/8266/2003", "A/Louisiana/4/2003", "A/Lousiana/4/2003", "A/OSAKA/31/2005",
        "A/Sari/388/2006", "A/HongKong/HK1/2008", "A/HongKong/HK1MA21-1/2008","A/HongKong/HK1MA21-4/2008", "A/HongKong/HK1MA21-2/2008",
        "A/HongKong/HK1MA21-3/2008", "A/HongKong/HK2/2008", "A/HongKong/HK2MA21-1/2008",
        "A/HongKong/HK2MA21-2/2008", "A/HongKong/HK2MA21-3/2008", "A/HongKong/HK4/2008",
        "A/HongKong/HK5/2008", "A/HongKong/HK5MA21-1/2008", "A/HongKong/HK5MA21-3/2008",
        "A/HongKong/HK6/2008", "A/HongKong/HK6MA21-2/2008", "A/HongKong/HK6MA21-3/2008",
        "A/HongKong/HK7/2008", "A/HongKong/HK8/2008", "A/HongKong/HK8MA21-1/2008",
        "A/HongKong/HK8MA21-2/2008", "A/HongKong/HK8MA21-3/2008", "A/HongKong/HK8MA21-4/2008",
        "A/HongKong/HK9/2008", "A/HongKong/HK9MA21-1/2008", "A/HongKong/HK9MA21-2/2008",
        "A/HongKong/HK9MA21-3/2008", "A/HongKong/HK10/2008", "A/HongKong/HK10MA21-1/2008",
        "A/HongKong/HK10MA21-2/2008", "A/HongKong/HK10MA21-3/2008", "A/HongKong/HK10MA21-4/2008",
        "A/HongKong/HK11MA21-1/2008", "A/HongKong/HK11MA21-3/2008", "A/HongKong/HK11MA21-4/2008",
        "A/HongKong/HK12/2008", "A/HongKong/HK12MA21-2/2008", "A/HongKong/HKMA12/2008",
        "A/HongKong/HKMA12A/2008", "A/HongKong/HKMA12B/2008", "A/HongKong/HKMA12D/2008",
        "A/HongKong/HKMA12E/2008", "A/HongKong/HKMA20B/2008", "A/HongKong/HKMA20E/2008", "A/Kansas/13/2009",
        "A/Busan/15453/2009", "A/Pennsylvania/14/2010", "A/Pennsylvania/40/2010", "A/Guam/AF2771/2011",
        "A/Indiana/8/2011", "A/Kenya/155/2011", "A/Kenya/168/2011", "A/Kenya/170/2011", "A/Nepal/142/2011",
        "A/Pennsylvania/09/2011", "A/Pennsylvania/9/2011", "A/Quebec/167936/2011", "A/Quebec/170658/2011",
        "A/India/6352/2012", "A/Indiana/08/2012", "A/Indiana/13/2012", "A/Ohio/34/2012",
        "A/Helsinki/942/2013", "A/Indiana/06/2013", "A/Indiana/11/2013", "A/Indiana/21/2013",
        "A/Jiangsu-Tianning/1707/2013", "A/HuNan/01/2014", "A/Jiangsu-Chongchuan/1830/2014",
        "A/Jiangsu-Chongchuan/12179/2014", "A/Ohio/2/2014", "A/Ohio/4319/2014", "A/SaoPaulo/3-34891/2014",
        "A/Wisconsin/24/2014", "A/NewJersey/53/2015", "A/SaoPaulo/36178/2015", "A/SaoPaulo/61282/2015",
        "A/SaoPaulo/65194/2015", "A/Michigan/39/2015", "A/Sydney/53/2015", "A/Michigan/82/2016",
        "A/Michigan/83/2016", "A/Michigan/84/2016", "A/Michigan/87/2016", "A/Michigan/89/2016",
        "A/Michigan/90/2016", "A/Michigan/91/2016", "A/Michigan/93/2016", "A/Michigan/94/2016",
        "A/Michigan/95/2016", "A/Michigan/96/2016", "A/Ohio/27/2016", "A/Ohio/28/2016", "A/Ohio/32/2016",
        "A/Ohio/33/2016", "A/Ohio/35/2016", "A/Zhejiang-Wuxin/1300/2016", "A/Nanjing/1/2010", "A/StPetersburg/5/2009",
        "A/Cambodia/NHRCC00001/2009","A/Cambodia/NHRCC00002/2009", "A/Cambodia/NHRCC00003/2009",
        "A/Cambodia/NHRCC00006/2009", ],
'h1n1pdm': ["A/Kenya/264/2012", "A/Iowa/39/2015", "A/Asturias/RR6898/2010", "A/Wisconsin/28/2011",
            "A/Brest/1161/2014", "A/Tomsk/273-MA1/2010", "A/Minnesota/46/2015", "A/Poland/16/2013",
            "A/Hungary/02/2013", "A/Hungary/16/2013", "A/California/07/2009NYMC-X18113/198",
            "A/Christchurch/16/2010NIB-74xp13/202", "A/Bari/166/2016", "A/Bari/167/2016", "A/Dakar/3/2014",
            "A/Arkansas/15/2013", "A/Wisconsin/87/2005", 'A/Norway/1929/2014', 'A/Ohio/9/2015'],
'vic':["B/Bangkok/SI17/2012", "B/Bangkok/SI58/2012", "B/Kol/2024/2008", "B/Kolkata/2024/2008"],
"yam":[]
}


clade_designations = {"h3n2":{
                           "3c3.a":[('HA1',128,'A'), ('HA1',142,'G'), ('HA1',159,'S')],
                           "3c3":   [('HA1',128,'A'), ('HA1',142,'G'), ('HA1',159,'F')],
                           "3c2.a": [('HA1',144,'S'), ('HA1',159,'Y'), ('HA1',225,'D'), ('HA1',311,'H'), ('HA2',160,'N')],
                           "171K": [('HA1',144,'S'), ('HA1',159,'Y'), ('HA1',171,'K'), ('HA1',225,'D'), ('HA1',311,'H'), ('HA2',77,'V'), ('HA2',155,'E'), ('HA2',160,'N')],
                           "3c2":   [('HA1',144,'N'), ('HA1',159,'F'), ('HA1',225,'N'), ('HA2',160,'N'), ('HA1',142,'R')],
                           "3c3.b": [('HA1',83,'R'), ('HA1',261,'Q'), ('HA1',62,'K'),  ('HA1',122,'D')]
                        },
                       "h1n1pdm":{
                            '2': [('HA1', 125, 'N'), ('HA1', 134 ,'A'), ('HA1', 183, 'S'), ('HA1', 31,'D'), ('HA1', 172,'N'), ('HA1', 186,'T')],
                            '3': [('HA1', 134 ,'T'), ('HA1', 183, 'P')],
                            '4': [('HA1', 125, 'D'), ('HA1', 134 ,'A'), ('HA1', 183, 'S')],
                            '5': [('HA1', 87, 'N'), ('HA1', 205, 'K'), ('HA1', 216, 'V'), ('HA1', 149, 'L')],
                            '6': [('HA1', 185,'T'),  ('HA1', 97, 'N'), ('HA1', 197, 'A')],
                            '6c':[('HA1', 234,'I'),  ('HA1', 97, 'N'), ('HA1', 197, 'A'), ('HA1', 283,'E')],
                            '6b':[('HA1', 163,'Q'),  ('HA1', 256, 'T'), ('HA1', 197, 'A'), ('HA1', 283,'E')],
                            '7': [('HA1', 143,'G'),  ('HA1', 97, 'D'), ('HA1', 197, 'T')],
                            '8': [('HA1', 186,'T'),  ('HA1', 272,'A')],
                            '6b.1':[('HA1', 163,'Q'),  ('HA1', 256, 'T'), ('HA1', 197, 'A'), ('HA1', 283, 'E'), ('SigPep', 13, 'T'), ('HA1', 84, 'N'), ('HA1', 162, 'N')],
                            '6b.2':[('HA1', 163,'Q'),  ('HA1', 256, 'T'), ('HA1', 197, 'A'), ('HA1', 283, 'E'), ('HA2', 164, 'G'), ('HA1', 152, 'T'), ('HA2', 174, 'E')]
                       },
                       "yam":{
                            '2':  [('HA1', 48,'K'), ('HA1', 108, 'A'), ('HA1', 150, 'S')],
                            '3':  [('HA1', 48,'R'), ('HA1', 108, 'P'), ('HA1', 150, 'I')],
                            '3a': [('HA1', 37,'A'), ('HA1', 298, 'E'), ('HA1', 48,'R'), ('HA1', 105, 'P'), ('HA1', 150, 'I')],
                            '172Q': [('HA1', 48,'R'), ('HA1', 108, 'P'), ('HA1', 150, 'I'), ('HA1', 116, 'K'), ('HA1', 172, 'Q')]
                       },
                       "vic":{
                            '1A': [('HA1', 75,'K'), ('HA1', 58, 'L'), ('HA1', 165, 'K')],
                            '1B': [('HA1', 75,'K'), ('HA1', 58, 'P'), ('HA1', 165, 'K')],
                            '117V': [('HA1', 75,'K'), ('HA1', 58, 'L'), ('HA1', 165, 'K'), ('HA1', 129, 'D'), ('HA1', 117, 'V')]
                        }
                     }

class flu_process(process):
    """process influenza virus sequences in mutliple steps to allow visualization in browser
        * filtering and parsing of sequences
        * alignment
        * tree building
        * frequency estimation of clades and mutations
        * export as json
    """

    def __init__(self, *args, **kwargs):
        super(flu_process, self).__init__(*args, **kwargs)
        print("Initializing flu_process")
        if 'titer' in kwargs and kwargs['titer']:
            hi_prefix = '_' + kwargs['titer']
        else:
            hi_prefix = '_hi'

        if 'segment' in kwargs:
            self.segment=kwargs['segment']
        else:
            self.segment='ha'

        if self.segment=='ha':
            self.HI_titer_fname = self.input_data_path   +hi_prefix +'_titers.tsv'
            self.HI_strains_fname = self.input_data_path +hi_prefix +'_strains.tsv'
        else:
            self.HI_titer_fname = self.input_data_path[:-3]   +hi_prefix +'_titers.tsv'
            self.HI_strains_fname = self.input_data_path[:-3] +hi_prefix +'_strains.tsv'

    def subsample(self, sampling_threshold, all_regions=False, **kwargs):
        nregions = len(regions)
        if not all_regions:
            region_threshold = int(np.ceil(1.0*sampling_threshold/nregions))
        self.sequence_count_total = defaultdict(int)
        self.sequence_count_region = defaultdict(int)
        for seq in self.seqs.all_seqs.values():
            self.sequence_count_total[(seq.attributes['date'].year,
                                  seq.attributes['date'].month)]+=1
            self.sequence_count_region[(seq.attributes['region'],
                                  seq.attributes['date'].year,
                                  seq.attributes['date'].month)]+=1

        if all_regions:
            def sampling_category(x):
                return (x.attributes['date'].year,
                        x.attributes['date'].month)
            def threshold_func(x):
                return sampling_threshold
        else:
            def sampling_category(x):
                return (x.attributes['region'],
                        x.attributes['date'].year,
                        x.attributes['date'].month)
            def threshold_func(x):
                if self.sequence_count_total[(x[1], x[2])]<sampling_threshold:
                    return sampling_threshold
                else:
                    region_counts = sorted([self.sequence_count_region[(r, x[1], x[2])] for r in regions])
                    if region_counts[0]>region_threshold:
                        return region_threshold
                    else:
                        left_to_fill = sampling_threshold - nregions*region_counts[0]
                        thres = region_counts[0]
                        for ri, rc in zip(range(len(regions)-1, 0, -1), region_counts[1:]):
                            if left_to_fill - ri*(rc-thres)>0:
                                left_to_fill-=ri*(rc-thres)
                                thres = rc
                            else:
                                thres += left_to_fill/ri
                                break
                        return max(1,int(thres))


        # load HI titer count to prioritize sequences
        HI_titer_count = {}
        with myopen(self.HI_strains_fname,'r') as ifile:
            for line in ifile:
                strain, count = line.strip().split()
                HI_titer_count[strain]=int(count)

        def sampling_priority(seq):
            sname = seq.attributes['strain']
            if sname in HI_titer_count:
                pr = HI_titer_count[sname]
            else:
                pr = 0
            return pr + len(seq.seq)*0.0001 - 0.01*np.sum([seq.seq.count(nuc) for nuc in 'NRWYMKSHBVD'])

        self.seqs.subsample(category = sampling_category,
                            threshold=threshold_func,
                            priority=sampling_priority, **kwargs)


    def subsample_priority_region(self, sampling_threshold, priority_region='north_america', fraction=0.5, **kwargs):
        self.sequence_count_total = defaultdict(int)
        self.sequence_count_region = defaultdict(int)
        for seq in self.seqs.all_seqs.values():
            self.sequence_count_total[(seq.attributes['date'].year,
                                  seq.attributes['date'].month)]+=1
            self.sequence_count_region[(seq.attributes['region'],
                                  seq.attributes['date'].year,
                                  seq.attributes['date'].month)]+=1

        def sampling_category(x):
            return (x.attributes['region'],
                    x.attributes['date'].year,
                    x.attributes['date'].month)


        def threshold_func(x):
            if x[0]==priority_region:
                return int(sampling_threshold*fraction)
            else:
                nregions = len(regions)-1
                total_threshold_world = sampling_threshold*(1-fraction)
                region_threshold = int(np.ceil(1.0*total_threshold_world/nregions))
                region_counts = sorted([self.sequence_count_region[(r, x[1], x[2])]
                                        for r in regions if r!=priority_region])
                if region_counts[0]>region_threshold:
                    return region_threshold
                else:
                    left_to_fill = total_threshold_world - nregions*region_counts[0]
                    thres = region_counts[0]
                    for ri, rc in zip(range(nregions-1, 0, -1), region_counts[1:]):
                        if left_to_fill - ri*(rc-thres)>0:
                            left_to_fill-=ri*(rc-thres)
                            thres = rc
                        else:
                            thres += left_to_fill/ri
                            break
                    return max(1,int(thres))


        # load HI titer count to prioritize sequences
        HI_titer_count = {}
        with myopen(self.HI_strains_fname,'r') as ifile:
            for line in ifile:
                strain, count = line.strip().split()
                HI_titer_count[strain]=int(count)

        def sampling_priority(seq):
            sname = seq.attributes['strain']
            if sname in HI_titer_count:
                pr = HI_titer_count[sname]
            else:
                pr = 0
            return pr + len(seq.seq)*0.0001 - 0.01*np.sum([seq.seq.count(nuc) for nuc in 'NRWYMKSHBVD'])

        self.seqs.subsample(category = sampling_category,
                            threshold=threshold_func,
                            priority=sampling_priority, **kwargs)


    def parse_age(self):
        for seq in self.seqs.all_seqs.values():
            tmp = seq.attributes['age']
            if tmp=='?':
                seq.attributes['age'] = np.nan
            elif tmp[-1]=='y':
                seq.attributes['age'] = float(tmp[:-1])
            elif tmp[-1]=='m':
                seq.attributes['age'] = float(tmp[:-1])/12.0
            elif tmp[-1]=='d':
                seq.attributes['age'] = float(tmp[:-1])/365.0
            else:
                try:
                    seq.attributes['age'] = float(tmp)
                except:
                    seq.attributes['age'] = np.nan


    def HI_model(self, **kwargs):
        '''
        estimate a tree and substitution model using titers titer_fname.
        '''
        from base.titer_model import tree_model, substitution_model
        ## TREE MODEL
        self.HI_tree = tree_model(self.tree.tree, titer_fname = self.HI_titer_fname, **kwargs)
        self.HI_tree.prepare(**kwargs)
        self.HI_tree.train(**kwargs)
        # add tree attributes to the list of attributes that are saved in intermediate files
        for n in self.tree.tree.find_clades():
            n.attr['cTiter'] = n.cTiter
            n.attr['dTiter'] = n.dTiter

        # SUBSTITUTION MODEL
        self.HI_subs = substitution_model(self.tree.tree, titer_fname = self.HI_titer_fname,**kwargs)
        self.HI_subs.prepare(**kwargs)
        self.HI_subs.train(**kwargs)


    def HI_export(self):
        from base.io_util import write_json
        prefix = self.build_data_path
        if hasattr(self, 'HI_tree'):
            # export the raw titers
            hi_data = self.HI_tree.compile_titers()
            write_json(hi_data, prefix+'titers.json')
            # export the tree model (avidities and potencies only)
            tree_model = {'potency':self.HI_tree.compile_potencies(),
                          'avidity':self.HI_tree.compile_virus_effects(),
                          'dTiter':{n.clade:n.dTiter for n in self.tree.tree.find_clades() if n.dTiter>1e-6}}
            write_json(tree_model, prefix+'titer_tree_model.json')
        else:
            print('Tree model not yet trained')

        if hasattr(self, 'HI_tree'):
            # export the substitution model
            subs_model = {'potency':self.HI_subs.compile_potencies(),
                          'avidity':self.HI_subs.compile_virus_effects(),
                          'substitution':self.HI_subs.compile_substitution_effects()}
            write_json(subs_model, prefix+'titer_subs_model.json')
        else:
            print('Substitution model not yet trained')


def H3N2_scores(tree, epitope_mask_version='wolf'):
    '''
    takes a H3N2 HA tree and assigns H3 specific characteristics to
    internal and external nodes
    '''
    def epitope_sites(aa):
        return aa[epitope_mask[:len(aa)]]

    def nonepitope_sites(aa):
        return aa[~epitope_mask[:len(aa)]]

    def receptor_binding_sites(aa):
        '''
        Receptor binding site mutations from Koel et al. 2014
        These are (145, 155, 156, 158, 159, 189, 193) in canonical HA numbering
        need to subtract one since python arrays start at 0
        '''
        sp = 16
        rbs = map(lambda x:x+sp-1, [145, 155, 156, 158, 159, 189, 193])
        return np.array([aa[pos] for pos in rbs])

    def get_total_peptide(node):
        '''
        the concatenation of signal peptide, HA1, HA1
        '''
        return np.fromstring(node.translations['SigPep']+node.translations['HA1']
                       + node.translations['HA2'], 'S1')

    def epitope_distance(aaA, aaB):
        """Return distance of sequences aaA and aaB by comparing epitope sites"""
        epA = epitope_sites(aaA)
        epB = epitope_sites(aaB)
        distance = np.sum(epA!=epB)
        return distance

    def nonepitope_distance(aaA, aaB):
        """Return distance of sequences aaA and aaB by comparing non-epitope sites"""
        neA = nonepitope_sites(aaA)
        neB = nonepitope_sites(aaB)
        distance = np.sum(neA!=neB)
        return distance

    def receptor_binding_distance(aaA, aaB):
        """Return distance of sequences aaA and aaB by comparing receptor binding sites"""
        neA = receptor_binding_sites(aaA)
        neB = receptor_binding_sites(aaB)
        distance = np.sum(neA!=neB)
        return distance

    epitope_map = {}
    with open('flu/metadata/h3n2_epitope_masks.tsv') as f:
        for line in f:
            (key, value) = line.strip().split()
            epitope_map[key] = value
    if epitope_mask_version in epitope_map:
        epitope_mask = np.fromstring(epitope_map[epitope_mask_version], 'S1')=='1'
    root = tree.root
    root_total_aa_seq = get_total_peptide(root)
    for node in tree.find_clades():
        total_aa_seq = get_total_peptide(node)
        node.attr['ep'] = epitope_distance(total_aa_seq, root_total_aa_seq)
        node.attr['ne'] = nonepitope_distance(total_aa_seq, root_total_aa_seq)
        node.attr['rb'] = receptor_binding_distance(total_aa_seq, root_total_aa_seq)


def make_date_ticks(ax, fs=12):
    from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
    years    = YearLocator()
    months = MonthLocator(range(1, 13), bymonthday=1, interval=2)
    yearsFmt = DateFormatter('%Y')
    monthsFmt = DateFormatter("%b")
    ax.tick_params(axis='x', which='major', labelsize=fs, pad=20)
    ax.tick_params(axis='x', which='minor', pad=7)
    ax.xaxis.set_major_locator(years)
    ax.xaxis.set_major_formatter(yearsFmt)
    ax.xaxis.set_minor_locator(months)
    ax.xaxis.set_minor_formatter(monthsFmt)


def plot_trace(ax, pivots, freq, err, n_std_dev=1, err_smoothing=3, show_errorbars=True, c='r', ls='-', label=None):
    ax.plot(pivots, freq, c=c, ls=ls, label=label)
    if show_errorbars:
        smerr = 1.0/np.convolve(1.0/err, np.ones(err_smoothing, dtype=float)/err_smoothing, mode='same')
        ax.fill_between(pivots, np.maximum(0,freq-n_std_dev*smerr),
                np.minimum(1,freq+n_std_dev*smerr),
                facecolor=c, linewidth=0, alpha=0.1)

def pivots_to_dates(pivots):
    date_bins = []
    for b in pivots:
        year = int(np.floor(b))
        month = int(round(12*(b%1)))
        if not month:
            month=12
            year-=1
        date_bins.append(datetime.strptime("%d-%d"%(year, month), "%Y-%m") - timedelta(days=8))
    return date_bins

def plot_frequencies(flu, gene, mutation=None, plot_regions=None, all_muts=False, ax=None, **kwargs):
    import seaborn as sns
    sns.set_style('whitegrid')
    cols = sns.color_palette()
    linestyles = ['-', '--', '-.', ':']
    if plot_regions is None:
        plot_regions=regions
    pivots = flu.pivots
    if ax is None:
        plt.figure()
        ax=plt.subplot(111)
    if type(mutation)==int:
        mutations = [x for x,freq in flu.mutation_frequencies[('global', gene)].iteritems()
                     if (x[0]==mutation)&(freq[0]<0.5 or all_muts)]
    elif mutation is not None:
        mutations = [mutation]
    else:
        mutations=None

    if mutations is None:
        for ri, region in enumerate(plot_regions):
            count=flu.mutation_frequency_counts[region]
            plt.plot(pivots, count, c=cols[ri%len(cols)], label=region)
    else:
        print("plotting mutations", mutations)
        for ri,region in enumerate(plot_regions):
            for mi,mut in enumerate(mutations):
                if mut in flu.mutation_frequencies[(region, gene)]:
                    freq = flu.mutation_frequencies[(region, gene)][mut]
                    err = flu.mutation_frequency_confidence[(region, gene)][mut]
                    c=cols[ri%len(cols)]
                    label_str = str(mut[0]+1)+mut[1]+', '+region
                    plot_trace(ax, pivots, freq, err, c=c,
                        ls=linestyles[mi%len(linestyles)],label=label_str, **kwargs)
                else:
                    print(mut, 'not found in region',region)
    ax.ticklabel_format(useOffset=False)
    ax.legend(loc=2)

def plot_sequence_count(flu, fname=None, fs=12):
    # make figure with region counts
    import seaborn as sns
    date_bins = pivots_to_dates(flu.pivots)
    sns.set_style('ticks')
    cols = sns.color_palette(n_colors=len(region_groups))
    region_label = {'global': 'Global', 'NA': 'N America', 'AS': 'Asia', 'EU': 'Europe', 'OC': 'Oceania'}
    fig, ax = plt.subplots(figsize=(8, 3))
    count_by_region = flu.mutation_frequency_counts
    drop = 3
    tmpcounts = np.zeros(len(flu.pivots[drop:]))
    plt.bar(date_bins[drop:], count_by_region['global'][drop:], width=18, \
            linewidth=0, label="Other", color="#bbbbbb", clip_on=False)
    for c,region in zip(cols, region_groups):
        if region!='global':
            plt.bar(date_bins[drop:], count_by_region[region][drop:],
                    bottom=tmpcounts, width=18, linewidth=0,
                    label=region_label[region], color=c, clip_on=False)
            tmpcounts += count_by_region[region][drop:]
    make_date_ticks(ax, fs=fs)
    ax.set_ylabel('Sample count')
    ax.legend(loc=3, ncol=1, bbox_to_anchor=(1.02, 0.53))
    plt.subplots_adjust(left=0.1, right=0.82, top=0.94, bottom=0.22)
    sns.despine()
    if fname is not None:
        plt.savefig(fname)



if __name__=="__main__":
    import argparse
    import matplotlib.pyplot as plt

    parser = argparse.ArgumentParser(description='Process virus sequences, build tree, and prepare of web visualization')
    parser.add_argument('-y', '--years_back', type = str, default = 3, help='number of years back to sample')
    parser.add_argument('--resolution', type = str, help ="outfile suffix, can determine -v and -y")
    parser.add_argument('-v', '--viruses_per_month_seq', type = int, default = 10, help='number of viruses sampled per month')
    parser.add_argument('-w', '--viruses_per_month_tree', type = int, default = 10, help='number of viruses sampled per month')
    parser.add_argument('-r', '--raxml_time_limit', type = float, default = 1.0, help='number of hours raxml is run')
    parser.add_argument('-d', '--download', action='store_true', default = False, help='load from database')
    parser.add_argument('-t', '--time_interval', nargs=2, help='specify time interval rather than use --years_back')
    parser.add_argument('-l', '--lineage', type = str, default = 'h3n2', help='flu lineage to process')
    parser.add_argument('--new_auspice', default = False, action="store_true", help='file name for new augur')
    parser.add_argument('--confidence', default = False, action="store_true", help='evaluate confidence intervals of internal node timing')
    parser.add_argument('--sampling', default = 'even', type=str,
                        help='sample evenly (even), or prioritize one region (region), otherwise sample randomly')
    parser.add_argument('--HI', default = 'hi', type=str, help='specify titer data to use')
    parser.add_argument('--load', action='store_true', help = 'recover from file')
    parser.add_argument('--no_tree', default=False, action='store_true', help = "don't build a tree")
    params = parser.parse_args()

    # default values for --viruses_per_month and --years_back from resolution
    if params.resolution == "2y":
        params.viruses_per_month_tree = 9*len(regions)
        params.viruses_per_month_seq = 10*len(regions)
        params.years_back = 2
    elif params.resolution == "3y":
        params.viruses_per_month_tree = 7*len(regions)
        params.viruses_per_month_seq = 10*len(regions)
        params.years_back = 3
    elif params.resolution == "6y":
        params.viruses_per_month_tree = 30
        params.viruses_per_month_seq = 100
        params.years_back = 6
    elif params.resolution == "12y":
        params.viruses_per_month_tree = 20
        params.viruses_per_month_seq = 100
        params.years_back = 12

    # construct time_interval from years_back
    if not params.time_interval:
        today_str = "{:%Y-%m-%d}".format(datetime.today())
        date_str = "{:%Y-%m-%d}".format(datetime.today() - timedelta(days=365.25 * params.years_back))
        params.time_interval = [date_str, today_str]

    if params.new_auspice:
        fname_prefix = "flu_"+params.lineage
    else:
        fname_prefix = params.lineage

    if params.sampling!="even":
        fname_prefix+='_'+params.sampling

    if params.HI!="hi":
        fname_prefix+='_'+params.HI

    input_data_path = '../fauna/data/'+params.lineage
    if params.resolution:
        store_data_path = 'store/'+ fname_prefix + '_' + params.resolution +'_'
        build_data_path = 'build/'+ fname_prefix + '_' + params.resolution +'_'
    else:
        store_data_path = 'store/'+ fname_prefix + '_'
        build_data_path = 'build/'+ fname_prefix + '_'

    ppy = 12
    flu = flu_process(input_data_path = input_data_path, store_data_path = store_data_path,
                   build_data_path = build_data_path, reference='flu/metadata/'+params.lineage+'_outgroup.gb',
                   proteins=['SigPep', 'HA1', 'HA2'], titer=params.HI,
                   method='SLSQP', inertia=np.exp(-1.0/ppy), stiffness=2.*ppy)


    if params.load:
        flu.load()
        flu.export()
    else:
        flu.load_sequences(fields={0:'strain', 2:'isolate_id', 3:'date', 4:'region',
                             5:'country', 7:"city", 8:"passage",9:'lab', 10:'age', 11:'gender'})

        time_interval = [datetime.strptime(x, '%Y-%m-%d').date() for x in params.time_interval]

        pivots = np.arange(time_interval[0].year+(time_interval[0].month-1)/12.0,
                           time_interval[1].year+time_interval[1].month/12.0, 1.0/ppy)
        flu.seqs.filter(lambda s:
            s.attributes['date']>=time_interval[0] and s.attributes['date']<time_interval[1])
        flu.seqs.filter(lambda s: len(s.seq)>=900)
        flu.seqs.filter(lambda s: s.name not in outliers[params.lineage])

        if params.sampling=='even':
            flu.subsample(params.viruses_per_month_seq, all_regions=False)
        elif params.sampling in regions:
            flu.subsample_priority_region(params.viruses_per_month_seq,
                                          priority_region=params.sampling, fraction=0.5)
        else:
            flu.subsample(params.viruses_per_month_seq, all_regions=True)

        flu.align()
        flu.dump()
        # first estimate frequencies globally, then region specific
        flu.estimate_mutation_frequencies(region="global", pivots=pivots, min_freq=.10)
        for region in region_groups.iteritems():
            flu.estimate_mutation_frequencies(region=region, min_freq=.10)

        if not params.no_tree:
            if params.sampling=='even':
                flu.subsample(params.viruses_per_month_seq, all_regions=False, repeated=True)
            elif params.sampling in regions:
                flu.subsample_priority_region(params.viruses_per_month_tree, priority_region=params.sampling,
                                              fraction=0.5, repeated=True)
            else:
                flu.subsample(params.viruses_per_month_seq, all_regions=True, repeated=True)

            flu.align()
            flu.build_tree()
            flu.clock_filter(n_iqd=3, plot=False, remove_deep_splits=True)
            flu.annotate_tree(Tc=0.03, timetree=True, reroot='best', confidence=params.confidence)
            for geo in geo_attributes:
                flu.tree.geo_inference(geo)

            flu.estimate_tree_frequencies()
            for region in regions:
                flu.estimate_tree_frequencies(region=region)
            flu.dump()

            flu.HI_model(criterium = lambda x:len(x.aa_mutations['HA1'])>0)
            H3N2_scores(flu.tree.tree)
            flu.dump()
            flu.matchClades(clade_designations[params.lineage])
            flu.export(extra_attr=['serum'], controls=attribute_nesting, geo_attributes=geo_attributes,
                       color_options=color_options, panels=panels)
            flu.HI_export()

    plot_sequence_count(flu, store_data_path+'tip_count.png')

    # plot an approximate skyline
    # from matplotlib import pyplot as plt
    # T = flu.tree.tt
    # plt.figure()
    # skyline, confidence = T.merger_model.skyline_inferred(gen = 50, confidence=2.0)
    # plt.fill_between(skyline.x, confidence[0], confidence[1], color=(0.8, 0.8, 0.8))
    # plt.plot(skyline.x, skyline.y)
    # plt.yscale('log')
    # plt.ylabel('effective population size')
    # plt.xlabel('year')
    # plt.ticklabel_format(axis='x',useOffset=False)
    # plt.savefig('%s_%s_skyline.png'%(params.lineage, params.resolution))

from __future__ import division, print_function
import sys
sys.path.append('')  # need to import from base
from base.io_util import make_dir, remove_dir, tree_to_json, write_json, myopen
from base.sequences import sequence_set, num_date
from base.tree import tree
from base.process import process
from datetime import datetime
import os
from glob import glob


#### Provide metadata ####
regions = ["japan_korea","china","south_asia","southeast_asia", "south_pacific","oceania","subsaharan_africa","west_asia","caribbean","south_america","central_america","north_america","europe"]
colors =  ["#4B26B1", "#3F4ACA", "#4272CE", "#4D92BF", "#5DA8A3", "#74B583", "#8EBC66", "#ACBD51", "#C8B944", "#DDA93C", "#E68B35", "#E3602D", "#DC2F24"]
region_cmap = zip(regions, colors)

color_options = {
    "region":{"key":"region", "legendTitle":"Region", "menuItem":"region", "type":"discrete", "color_map": region_cmap},
    "num_date":{"key":"num_date", "legendTitle":"Sampling date", "menuItem":"date", "type":"continuous"},
}

attribute_nesting = {'geographic location':['region'], 'authors':['authors']}
panels = ['tree', 'map', 'entropy']#, 'frequencies']

##### Utils #####
def select_serotype(infile, path, serotype):
    '''
    Takes all-serotype fasta file, path to save output, and desired serotype.
    Writes appropriate subset of sequences to dengue_serotype.fasta
    Returns path to output file as string.
    '''
    from Bio import SeqIO
    sequences = [ i for i in SeqIO.parse(infile, 'fasta') if i.description.split('/')[0] == 'DENV%s'%serotype ]
    SeqIO.write(sequences, path+'dengue_%s.fasta'%serotype, 'fasta')
    return path+'dengue_%s.fasta'%serotype

##### The 'heavy lifting' #####
class dengue_process(process):
    def __init__(self, **kwargs):
        '''
        Congruent with other nextstrain builds, dengue_process is a catch-all class
        that initially holds the input data paths and params arguments.
        '''
        super(process, self).__init__()

        ##### Handle serotype-specific file input/output. #####
        self.serotype = kwargs['serotype']
        self.lineage = 'dengue_%s'%self.serotype
        if self.serotype == 'all': # For all-serotype build, use dengue 4 outgroup and look for files like dengue.fasta
            self.reference_fname = './dengue/metadata/dengue_denv4_outgroup.gb'
            newest_sequence_file = sorted(glob('../fauna/data/%s.fasta'%self.lineage), key=lambda f: os.path.getmtime(f))[-1]
        else:
            self.reference_fname = './dengue/metadata/%s_outgroup.gb'%self.lineage
            try: # Look for a serotype-specific fasta
                newest_sequence_file = sorted(glob('../fauna/data/%s*.fasta'%self.lineage), key=lambda f: os.path.getmtime(f))[-1]
            except: # If it doesn't exist, try to pull serotype-specific sequences out of the all-serotype fasta (warn the user of this behavior)
                newest_sequence_file = select_serotype('../fauna/data/dengue_all.fasta', '../fauna/data/', self.serotype)
                print('WARNING: Did not find serotype-specific fasta file.\nPulled sequences with serotype %s from all-serotype fasta file %s\nWrote these to file %s'%(self.serotype, '../fauna/data/dengue.fasta', newest_sequence_file))

        self.input_data_path = newest_sequence_file.split('.fasta')[0]
        self.sequence_fname = newest_sequence_file
        self.store_data_path = 'store/'+self.lineage + '_'
        self.build_data_path = 'build/'+self.lineage + '_'
        self.proteins = ['C', 'M', 'E', 'NS1', 'NS2A', 'NS2B', 'NS3', 'NS4A', '2K', 'NS4B', 'NS5']

        ##### Initialize process object #####
        self.dengue = process(input_data_path = self.input_data_path,
                       store_data_path = self.store_data_path,
                       build_data_path = self.build_data_path,
                       proteins=self.proteins,
                       reference=self.reference_fname,
                       method='SLSQP',
                       lat_long_fname='../fauna/source-data/geo_lat_long.tsv')

    def run_build(self, steps):
            '''
            Process then also gets the alignment, tree, and entropy objects.
            '''

            if 'load' in steps: ### Load checkpoint from file to bypass subsampling/alignments/tree building for time.
                print('\nLoading from file.....\n\n')
                self.dengue.load()

            if 'align' in steps: ### Read in data, align sequences, deal with metadata.
                print('\nSubsampling, filtering, and aligning sequences.....\n\n')

                self.fasta_fields = {0:'strain', 1:'accession', 2:'date', 3:'region', 4:'country',
                                5:'division', 6: 'location', 7: 'authors', 8: 'url'}
                self.dengue.load_sequences(fields=self.fasta_fields)
                assert self.dengue.seqs != None, 'ERROR: No sequences'
                self.dropped_strains = ['DENV1/VIETNAM/BIDV992/2006', 'DENV1/FRANCE/00475/2008', 'DENV1/VIETNAM/BIDV3990/2008', 'DENV2/HAITI/DENGUEVIRUS2HOMOSAPIENS1/2016', # Probable recombinants
                                        'DENV2/AUSTRALIA/QML22/2015'] # Suspiciously far diverged
                self.dengue.seqs.filter(lambda s: len(s.seq)>=5000, leave_ref=True) # Toss sequences shorter than 5000 bases
                self.dengue.seqs.filter(lambda s: s.id not in self.dropped_strains, leave_ref=True) # Toss sequences in dropped_strains
                self.dengue.seqs.filter(lambda s: s.attributes['region'] not in ['', '?'], leave_ref=True) # Toss sequences without location data
                self.dengue.seqs.subsample(category = lambda x:(x.attributes['region'], # Subsample per region per month
                                                         x.attributes['date'].year,
                                                         x.attributes['date'].month), threshold=params.viruses_per_month)
                for s in self.dengue.seqs.seqs.values():
                    s.attributes['url'] = 'https://www.ncbi.nlm.nih.gov/nuccore/%s'%s.attributes['accession'] # Add genbank URL to each sequence
                self.dengue.align(debug=True) # Align with mafft; save the .fasta alignment
                self.dengue.dump() # Save a checkpoint

            if 'tree' in steps: ### Build a tree (fast tree --> raxml)
                print('\nBuilding a tree.....\n\n')
                assert self.dengue.seqs != None, 'ERROR: No sequences, cannot build tree.'
                self.dengue.build_tree(debug=True) # Save the newick from raxml
                self.dengue.dump() # Save a checkpoint


            if 'treetime' in steps:
                print('\nStarting treetime clock filtering, rerooting, and tree annotation.....\n\n')
                assert self.dengue.seqs != None, 'ERROR: No sequences, cannot run treetime.'
                assert self.dengue.tree, 'ERROR: No tree, cannot run treetime.'
                self.dengue.clock_filter(n_iqd=3, plot=True) # Toss sequences that don't show a linear relationship between sampling time and root-to-tip divergence.
                self.dengue.annotate_tree(Tc=False, timetree=True, reroot='best') # Map attributes from the `fasta_headers` to the tree

            if 'geo' in steps: # Use a "mugration model" to trace how viruses have moved between regions
                print('\nRunning geo inference.....\n\n')
                assert 'treetime' in steps, 'ERROR: Must run step "treetime" to do geo inference.'
                self.dengue.tree.geo_inference('region')

            if 'export' in steps:
                print('\nExporting.....\n\n')
                assert self.dengue.tree, 'ERROR: No tree object to export'
                self.date_range = {'date_min': self.dengue.tree.getDateMin(), 'date_max': self.dengue.tree.getDateMax()} # Set default date range according to the tree height
                self.dengue.export(controls = attribute_nesting, geo_attributes = 'region', date_range = self.date_range, color_options=color_options, panels=panels, defaults={'geoResolution': 'region'}) # Export to JSON

##### Config #####
if __name__=="__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Process virus sequences, build tree, and prepare web visualization',
    epilog="Which build steps to run can be configured in __main__; valid options are: ['load', 'align', 'tree', 'treetime', 'geo', 'export']")
    parser.add_argument('-v', '--viruses_per_month', type = int, default = 3, help='Number of viruses sampled per region per month')
    parser.add_argument('-s', '--serotype', type = str.lower, choices=['denv1', 'denv2', 'denv3', 'denv4', 'all'], default='all', help = 'Which serotype of dengue to build trees for; all = include all serotypes in one build.')
    parser.add_argument('--load', action='store_true', help = 'Recover alignment and tree from file, then run treetime, geo, and export steps.')
    parser.add_argument('--run_multiple', action='store_true', help = 'Run all 5 builds (4 serotype-specific + 1 all-serotype)')
    params = parser.parse_args()

    if params.load:
        steps = ['load', 'treetime', 'geo', 'export']
    else:
        steps = ['align', 'tree', 'treetime', 'geo', 'export']

    if params.run_multiple or not params.serotype: # Run all 5 builds.
        print('\n\nRunning all 5 builds.\n\n')
        for s in ['denv1', 'denv2', 'denv3', 'denv4', 'all']:
            params.serotype=s
            output = dengue_process(**params.__dict__)
            print('\nInitializing these steps for serotype %s:\n'%s, steps, '\n\n')
            output.run_build(steps)
    else: # Run single build and quit.
        output = dengue_process(**params.__dict__)
        print('\nInitializing these steps for serotype %s:\n'%params.serotype, steps, '\n\n')
        output.run_build(steps)

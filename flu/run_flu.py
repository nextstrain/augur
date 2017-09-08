import os

prepare = False

for lineage in ['h3n2', 'h1n1pdm', 'vic', 'yam']:
	for resolution in ['2y', '3y', '6y']:
		for passage in ['cell', 'egg']:
			for titer in ['hi', 'fra']:
				if lineage!='h3n2' and titer=='fra':
					continue

				call = ['python', 'flu.prepare.py', '-r', resolution,  '-l', lineage, '--titers', '../../fauna/data/%s_cdc_%s_%s_titers.tsv'%(lineage, titer, passage),
					'--sequences', '../../fauna/data/%s.fasta'%lineage, '--file_prefix', 'flu_%s_ha_%s_%s_%s'%(lineage, resolution, passage, titer)]
				if passage=='cell' and (titer=='hi' or lineage!='h3n2'):
					call.append('--complete_frequencies')

				os.system(' '.join(call))
				call = ['qsub', 'submit_script.sh', 'flu.process.py', '-j', 'prepared/flu_%s_ha_%s_%s_%s.json'%(lineage, resolution, passage, titer)]
				if not (passage=='cell' and (titer=='hi' or lineage!='h3n2')):
					call.append('--no_mut_freqs')

				os.system(' '.join(call))
				print(' '.join(call))

import os

prepare = False

for lineage in ['h3n2', 'h1n1pdm', 'vic', 'yam']:
	for resolution in ['2y', '3y', '6y']:
		for passage in ['cell', 'egg']:
			for titer in ['hi', 'fra']:
				if lineage!='h3n2' and titer=='fra':
					continue

				if prepare:
					call = ['python', 'flu.prepare.py', '-r', resolution,  '-l', lineage, '--titers', '../../fauna/data/%s_cdc_%s_%s_titers.tsv'%(lineage, titer, passage),
						'--sequences', '../../fauna/data/%s.fasta'%lineage, '--dataset', '%s_%s'%(titer, passage)]
				else:
					call = ['qsub', 'submit_script.sh', 'flu.process.py', '-j', 'prepared/flu_%s_ha_%s_%s_%s.json'%(lineage, resolution, titer, passage)]

				print(' '.join(call))
				os.system(' '.join(call))


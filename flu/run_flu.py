import os


for lineage in ['h3n2', 'h1n1pdm', 'vic', 'yam']:
	for resolution in ['2y', '3y', '6y']:
		for passage in ['cell', 'egg']:
			for titer in ['hi', 'fra']:
				if lineage!='h3n2' and titer=='fra':
					continue

				call = ['python', 'flu.prepare.py', '-l', lineage, '--titers', '../../fauna/data/%s_cdc_%s_%s_titers.tsv'%(lineage, titer, passage),
						'--sequences', '../../fauna/data/%d.fasta'%lineage, '--dataset', '%s_%s'%(titer, passage)]

				os.system(' '.join(call))

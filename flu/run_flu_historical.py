import os

for lineage in ['h3n2', 'h1n1pdm', 'vic', 'yam']:
    for time_interval in [('%d-01-01'%y, '%d-12-31'%(y+9)) for y in [1995, 2000,2005, 2010]]:
        call = ['python', 'flu.prepare.py', '--time_interval',time_interval[0], time_interval[1],
                '-l', lineage, '--titers', '../../fauna/data/%s_crick_hi_cell_titers.tsv'%(lineage),
                '--sequences', '../../fauna/data/%s.fasta'%lineage]
        print(' '.join(call))
        os.system(' '.join(call))

        call = ['qsub', 'submit_script.sh', 'flu.process.py','--pivot_spacing 3', '-j',
                'prepared/flu_%s_ha_%s_%s.json'%(lineage, time_interval[0], time_interval[1])]
        print(' '.join(call))
        os.system(' '.join(call))


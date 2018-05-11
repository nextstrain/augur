import os

def build(
    system = "local",
    segments = ['pb2', 'pb1', 'pa', 'ha', 'np', 'na', 'mp', 'ns']
    ):
    call = ['python',
        'avian.prepare.py',
        '--segments', ' '.join(segments)
        ]
    print(' '.join(call))
    os.system(' '.join(call))
    for segment in segments:
        call = [
            'avian.process.py',
            '--json', 'prepared/flu_avian_h7n9_%s.json'%(segment)]
        if (system == "qsub"):
            call = ['qsub', 'submit_script.sh'] + call
        elif (system == "sbatch"):
            call = ['python'] + call
            concat = '"' + ' '.join(call) + '"'
            call = ['sbatch', '-n', '1', '-c', '2', '--mem', '16192', '--time', '12:00:00', '--wrap', concat]
        elif (system == "local"):
            call = ['python'] + call
        print(' '.join(call))
        os.system(' '.join(call))

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Run flu builds')
    parser.add_argument('-s', '--system', type = str, default = 'local', help='where to run, local, qsub or sbatch')
    parser.add_argument('--segments', nargs='+', type = str,  help ="segments to include")
    params = parser.parse_args()

    if params.segments is None:
        params.segments = ['pb2', 'pb1', 'pa', 'ha', 'np', 'na', 'mp', 'ns']

    build(system = params.system, segments = params.segments)

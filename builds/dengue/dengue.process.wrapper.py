import subprocess
import sys

serotypes = ['denv1', 'denv2', 'denv3', 'denv4', 'all']

if len(sys.argv) > 1:
	cmds = ["python ./dengue.process.py -s %s %s"%(s, ' '.join(sys.argv[1:])) for s in serotypes]
else:
	cmds = ["python ./dengue.process.py -s %s"%s for s in serotypes]

running = []
job_id = 0

while True:
	running = [r for r in running if r.poll() is None]
	if len(running) != 2:
		if job_id == len(cmds):
			break

		while len(running) < 2:
			proc = subprocess.Popen(cmds[job_id].split())
			print 'STARTING JOB ID %d'%job_id
			job_id += 1
			running.append(proc)

import os
import datetime
import sys
import inspect
from qualikiz.qualikizrun import Run, EmptyJob, run_job, generate_input
from qualikiz.edisonbatch import Srun, Batch
""" Creates a mini-job based on the reference example """

runsdir = ''
if len(sys.argv) == 2:
    if sys.argv[1] != '':
        runsdir = os.path.abspath(sys.argv[1])

rootdir =  os.path.abspath(
   os.path.join(os.path.abspath(inspect.getfile(inspect.currentframe())),
                '../../'))
run = Run(rootdir, runsdir=runsdir)
srun = Srun('QuaLiKiz', 24)
batch = Batch('debug', srun)
batch.name = 'mini'
batch.optimize_sbatch()
batch.maxtime = datetime.timedelta(seconds=1*60)
job = EmptyJob(run.rootdir, 'mini', 'QuaLiKiz', batch)
run.jobs.update({job.name: job})
run.to_file()

resp = input('Submit job to queue? [Y/n]')
if resp == '' or resp == 'Y' or resp == 'y':
    job_path = os.path.join(run.runsdir, batch.name)
    generate_input(job_path)
    run_job(job_path)

import datetime
import math
import sys
import os
import inspect
from warnings import warn

import numpy as np

from qualikiz.qualikizrun import Run, EmptyJob, run_jobs, generate_inputs
from qualikiz.edisonbatch import Srun, Batch
from qualikiz.inputfiles import Electron, Ion, IonList, QuaLiKizRun

""" Creates a scan over number of nodes and xpoints"""
runsdir = ''
if len(sys.argv) == 2:
    if sys.argv[1] != '':
        runsdir = os.path.abspath(sys.argv[1])

rootdir =  os.path.abspath(
   os.path.join(os.path.abspath(inspect.getfile(inspect.currentframe())),
                '../../'))
run = Run(rootdir, runsdir=runsdir)

# Load default file
with open(os.path.join(run.qualikizdir,
                       'parameters_template.py')) as file_:
    parameters_template_script = file_.readlines()

xpoints_scan = [18*48, 2*18*48, 3*18*48, 4*18*48, 6*18*48, 8*18*48]
prev_nice_potxpt = [i for i in range(1, xpoints_scan[-1])]
for xpoints in xpoints_scan:
   nice_potxpt = []
   for potxpt in range(1, 3*xpoints):
      if potxpt in prev_nice_potxpt:
         if (xpoints % potxpt) == 0:
            if (3*xpoints/potxpt % 24) == 0:
               nice_potxpt.append( potxpt)
   prev_nice_potxpt = nice_potxpt
   #print (nice_potxpt)

xpoints_per_tasks_scan = [9, 12, 18, 27, 36, 54]
xpoints_per_task_grid, xpoints_grid = np.meshgrid(xpoints_per_tasks_scan, xpoints_scan)
tasks_grid = (3 * xpoints_grid / xpoints_per_task_grid).astype(int)
print ('xpoints_per_task = ' + str(xpoints_per_tasks_scan))
print ('xpoints = ' + str(xpoints_scan))

nodes_grid = tasks_grid / 24
maxtime = datetime.timedelta(seconds=25*60)
max_cputime_grid = nodes_grid * 24 * maxtime.total_seconds() / 3600
worst_case_cputime = []


mask = max_cputime_grid > 100
mask = np.ones_like(max_cputime_grid)
mask[:2,:] = 0
mask[1,0:2] = 1

xpoints_per_task_grid = np.ma.array(xpoints_per_task_grid, mask=mask)
tasks_grid = np.ma.array(tasks_grid, mask=mask)
xpoints_grid = np.ma.array(xpoints_grid, mask=mask)
max_cputime_grid = np.ma.array(max_cputime_grid)
max_cputime_grid.__setmask__(mask)

print ('xpoints per task')
print (xpoints_per_task_grid)
print ('xpoints')
print (xpoints_grid)
print ('tasks')
print (tasks_grid)
print ('max CPUh')
print (max_cputime_grid)


print ('scanning ' + str(max_cputime_grid.count()) + ' sbatch parameters')

print ('uses at max ' + str(max_cputime_grid.sum()) + ' CPUh')



for tasks_slice, xpoints_slice in zip(tasks_grid, xpoints_grid):
    for tasks, xpoints in zip(tasks_slice.compressed(), xpoints_slice.compressed()):
        binary_basename = 'QuaLiKiz'
        parameters_script = []
        for line in parameters_template_script:
            if 'xpoints =' in line:
                parameters_script.append('xpoints = ' + str(xpoints) + '\n')
            elif 'npoints =' in line:
                npoints = 8
                parameters_script.append('npoints = ' + str(npoints) + '\n')
            #elif 'setup' in line:
            #    parameters_script.append('hypercube = run.setup_hypercube()\n')
            #elif 'to_file' in line:
            #    parameters_script.append('run.to_file(hypercube)\n')
            else:
                parameters_script.append(line)


        name = 'dimxn' + str(xpoints*3*npoints*2) + \
                'mpi' + str(tasks)
    
        if xpoints*3 % tasks != 0:
             warn('xpoints not divisable')
        srun = Srun(binary_basename, tasks)
        batch = Batch('regular', srun, HT=True)
        batch.name = 'OMP_MPI_' + name
        try:
            batch.optimize_sbatch()
        # pylint: disable=broad-except
        except Exception as exception:
            warn('Could not optimize batch settings')
            warn(exception)
        else:
            batch.maxtime = maxtime
            worst_case_cputime.append(batch.cores_per_node * batch.nodes * maxtime)
            job = EmptyJob(run.rootdir, name, binary_basename, batch, parameters_script=parameters_script)
            run.jobs.update({job.name: job})
run.to_file(overwrite=True)

print ('Generating input files', end='', flush=True)
generate_inputs(run.runsdir, dotprint=True)
print ('\n')
resp = input('Submit all jobs in rundir to queue? [Y/n]')
if resp == '' or resp == 'Y' or resp == 'y':
    run_jobs(run.runsdir)

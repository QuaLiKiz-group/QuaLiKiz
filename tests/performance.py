import datetime
import math
import sys
import os
import inspect
import copy
from warnings import warn
from collections import OrderedDict

import numpy as np

from qualikiz.qualikizrun import QuaLiKizBatch, QuaLiKizRun
from qualikiz.edisonbatch import Srun, Sbatch
from qualikiz.inputfiles import Electron, Ion, IonList, QuaLiKizPlan

""" Creates a scan over number of nodes and xpoints"""
dirname = 'runs'
if len(sys.argv) == 2:
    if sys.argv[1] != '':
        dirname = sys.argv[1]

# We know where this script lives, so we can find the rootdir like this
rootdir =  os.path.abspath(
   os.path.join(os.path.abspath(inspect.getfile(inspect.currentframe())),
                '../../'))

# Find some of the paths we need
runsdir = os.path.join(rootdir, dirname)
bindir = os.path.join(rootdir, 'QuaLiKiz+pat')
templatepath = os.path.join(rootdir, 'tools/qualikiz/parameters_template.json')

# Load the default QuaLiKiz plan we will use as base
qualikiz_plan_base = QuaLiKizPlan.from_json(templatepath)
print (qualikiz_plan_base['xpoint_base']['special']['kthetarhos'])
dimn = len(qualikiz_plan_base['xpoint_base']['special']['kthetarhos'])

# Set the amount of points in one scan. At dimxn=552960 we run out of memory
points_scan = [3*9*48, 2700, 15*9*48] 
xpoints_scan = [len(qualikiz_plan_base['scan_dict']) * points for points in points_scan]
print ('dimxn = ' + str(dimn * np.array(xpoints_scan)))

# Some legacy to find nicely splittable points
#dimensions = 4
#prev_nice_potxpt = [i for i in range(1, xpoints_scan[-1])]
#for xpoints in xpoints_scan:
#   nice_potxpt = []
#   for potxpt in range(1, dimensions*xpoints):
#      if potxpt in prev_nice_potxpt:
#         if (xpoints % potxpt) == 0:
#            if (dimensions*xpoints/potxpt % 24) == 0:
#               nice_potxpt.append( potxpt)
#   prev_nice_potxpt = nice_potxpt
#   print (nice_potxpt)

numnodes_scan = [32, 48, 60]
numcores_scan = [24 * numnodes for numnodes in numnodes_scan]
print (numcores_scan)
cores_grid, xpoints_grid = np.meshgrid(numcores_scan, xpoints_scan)

mask = np.zeros_like(cores_grid)
#mask[1:, 0] = 1
#mask[:1, 1] = 1

xpoints_grid = np.ma.array(xpoints_grid, mask=mask)
cores_grid = np.ma.array(cores_grid, mask=mask)

print ('xpoints')
print (xpoints_grid)
print ('cores')
print (cores_grid)
print ('xpoints per core')
print (np.array(xpoints_grid) / np.array(cores_grid))

batch_list = []
est_cpu_time = []
est_wall_time = []

for cores_slice, xpoints_slice in zip(cores_grid, xpoints_grid):
    for cores, xpoints in zip(cores_slice, xpoints_slice):
        if not np.ma.is_masked(cores) and not np.ma.is_masked(xpoints):
            qualikiz_plan = copy.deepcopy(qualikiz_plan_base)

            name = 'dimxn' + str(xpoints*dimn) + \
                    'cores' + str(cores)


            points = xpoints / len(qualikiz_plan['scan_dict'])
            new_scan_dict = OrderedDict()
            for key, value in qualikiz_plan['scan_dict'].items():
                new_scan_dict[key] = np.linspace(value[0], value[-1], points).tolist()

            qualikiz_plan['scan_dict'] = new_scan_dict
               
            binreldir = os.path.relpath(bindir,
                                        os.path.join(runsdir, name))
            run = QuaLiKizRun(runsdir, name,
                              binreldir,
                              qualikiz_plan=qualikiz_plan)
            est_cpu_time.append(run.estimate_cputime(cores))
            est_wall_time.append(run.estimate_walltime(cores))
            runlist = [run]
            try:
               batch = QuaLiKizBatch(runsdir, name, runlist, cores, partition='debug')
            except:
               warn('careful! Too long!')
            batch.prepare(overwrite_batch=True)
            batch_list.append(batch)
        else:
            est_cpu_time.append(0)
    
print ('CPUh')
print (np.reshape(est_cpu_time, cores_grid.shape)/3600)
print ('takes an estimated ' + str(np.sum(est_cpu_time)/3600) + ' CPUh')
print ('Wallh')
print (np.reshape(est_wall_time, cores_grid.shape)/60)
print ('takes an estimated ' + str(np.sum(est_wall_time)/60) + ' Wallm')

print ('Generating input files', end='', flush=True)
for batch in batch_list:
    batch.generate_input(dotprint=True)
print ('\n')
resp = input('Submit all jobs in runsdir to queue? [Y/n]')
if resp == '' or resp == 'Y' or resp == 'y':
   for batch in batch_list:
      batch.queue_batch()

""" Make a plot of the performance using a polldump """
import sys
import pickle
import os
from collections import OrderedDict
import re

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np

from qualikiz.tabulate import tabulate
from qualikiz.qualikizrun import acctstr_to_timedelta

if len(sys.argv) < 2:
    raise Exception('Please supply a file to analyse')
file_ = os.path.abspath(sys.argv[1])
with open(file_, 'rb') as file__:
    acctlist = pickle.load(file__)
results = []
# Extract the information we want
for entry in acctlist:
    result = OrderedDict()
    for job in entry:
        if job['State'] == 'COMPLETED':
            if job['JobID'].endswith('.0'):
                ntasks = int(job['NTasks'])
                result['tasks'] = ntasks
            elif job['JobID'].endswith('.batch'):
                pass
            else:
                name = job['JobName']
                # As we intelligently named our scripts, we can extract dimxn
                dimxn, __ = [int(n) for n in re.split('\D', name) if n is not '']
                result['dimxn'] = dimxn
                ncpus = int(job['AllocCPUS'])
                result['JobID'] = job['JobID']
                walltime = acctstr_to_timedelta(job['Elapsed'])
                NNodes = int(job['NNodes'])
                result['NNodes'] = NNodes
                rawCores = 24 * NNodes
                result['ncpus'] = ncpus
                result['rawcores'] = rawCores
                result['dimxnprawcores'] = dimxn / rawCores
                rawCPUs = walltime.seconds * rawCores
                rawCPUh = rawCPUs / 60 / 60
                if rawCPUs > 0:
                    result['rawCPUh'] = rawCPUh
                else:
                    result['rawCPUh'] = np.nan

                rawWalls = walltime.seconds
                rawWallh = rawWalls / 60 / 60
                if rawWalls > 0:
                    result['rawWallh'] = rawWallh
                else:
                    result['rawWallh'] = np.nan

                results.append(result)

results_sorted = OrderedDict()
for result in results:
    try:
        results_sorted[result['dimxn']].append(result)
    except KeyError:
        results_sorted[result['dimxn']] = [result]

for dimxn, result_batch in results_sorted.items():
    result_batch = sorted(result_batch, key=lambda t: t['NNodes'])
    print (tabulate.tabulate(result_batch, headers='keys') + '\n')


xs = np.empty([0])
ys = np.empty([0])
xind = 'dimxnprawcores'
yind = 'dimxn'
zind = 'rawWallh'
xs = np.empty([0], dtype=type(results[0][xind]))
ys = np.empty([0], dtype=type(results[0][yind]))
for result in results:
    if result[xind] not in xs:
        xs = np.append(xs, result[xind])
    if result[yind] not in ys:
        ys = np.append(ys, result[yind])
xs = np.sort(xs)
ys = np.sort(ys)

zs = np.full([len(xs), len(ys)], np.nan)
for result in results:
    ix = np.flatnonzero(xs == result[xind])
    iy = np.flatnonzero(ys == result[yind])
    zs[ix, iy] = result[zind]

ygrid, xgrid = np.meshgrid(ys, xs)
print (xgrid)
print (ygrid)
print (zs)

fig = plt.figure()
ax1 = fig.add_subplot(111)

for x_slice, y_slice, z_slice in zip(xgrid.T, ygrid.T, zs.T):
    ax1.loglog(x_slice[np.isfinite(z_slice)], z_slice[np.isfinite(z_slice)], basex=2, label=y_slice)
plt.legend(ys)
plt.legend([str(y) + ' points' for y in ys])
plt.gca().invert_xaxis()

ax1.set_xlabel('points per core')
ax1.set_ylabel('Walltime [h]')

fig = plt.figure()
ax2 = fig.add_subplot(111)
for x_slice, y_slice, z_slice in zip(xgrid, ygrid, zs):
    ax2.loglog(y_slice[np.isfinite(z_slice)], z_slice[np.isfinite(z_slice)], basex=2, label=x_slice)
plt.legend([str(x) + ' points per core' for x in xs])
ax2.set_xlabel('points')
ax2.set_ylabel('Walltime [h]')

plt.tight_layout()





plt.show()

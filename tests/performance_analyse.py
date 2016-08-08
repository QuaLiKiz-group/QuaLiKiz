""" Make a plot of the performance using a polldump """
import sys
import os
from warnings import warn
from itertools import islice
import sqlite3

import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np

from qualikiz.tabulate.tabulate import tabulate
from qualikiz.sacct import acctstr_to_timedelta

if len(sys.argv) < 2:
    raise Exception('Please supply a file to analyse')
file_ = os.path.abspath(sys.argv[1])

db = sqlite3.connect(file_)
db.execute('DROP TABLE IF EXISTS profile')
db.execute('''CREATE TABLE profile (
    Jobnumber  INTEGER,
    Dimxn      INTEGER,
    NumrawCPU  INTEGER,
    Pointspcore REAL
    )''')

# First calculate the parameters of interest
query = db.execute('''SELECT stdout.Jobnumber, Dimx, Dimn, Nodes, Tasks
                   FROM jobdata
                       JOIN stdout ON stdout.Jobnumber=jobdata.Jobnumber''')
for line in query:
    jobnumber = line[0]
    dimxn = line[1] * line[2]
    numrawCPU = line[3] * 24
    pointspercore = dimxn / numrawCPU

    db.execute('''
    INSERT INTO profile(
    Jobnumber, Dimxn, NumrawCPU, Pointspcore)
    VALUES (?, ?, ?, ?)''', (jobnumber, dimxn, numrawCPU, pointspercore))
db.commit()

# Sanity check, the total time should be the sum of the seperate parts
query = db.execute('''SELECT Jobnumber, (Total), First_MPI_AllReduce, (Eigenmodes), 
                      (Saturation), Second_MPI_AllReduce, (Initialization), (Output)
                   FROM stdout''')
for line in query:
    sanity = np.isclose(line[1], np.sum(line[1+1:]), atol=1e-3, rtol=1e-2)
    if not sanity:
        warn('Job ' + str(line[0]) + ' is insane! Sum of QuaLiKiz parts time != total time')
labels = ['First MPI AllReduce', 'Eigenmodes', 'Saturation', 'Second MPI Allreduce', 'Init', 'Output']

# If sacct table exists, also check reported time with NERSC system
query = db.execute("SELECT name FROM sqlite_master WHERE type='table'")
rows = query.fetchall()
tables = [x[0] for x in rows]
if 'sacct' in tables:
    query = db.execute('''SELECT JobID, Elapsed, NNodes, Total
                       FROM sacct
                           JOIN stdout ON sacct.JobID=stdout.Jobnumber''')
    headers = ['JobID', 'NERSC wall', 'NERSC nodes', 'NERSC CPU', 'QuaLiKiz wall']
    table = []
    tot_cpu = 0
    for line in query:
        calc_cpu = acctstr_to_timedelta(line[1]) * line[2] * 24
        tot_cpu += calc_cpu.total_seconds()
        table.append([line[0],
                      line[1],
                      line[2],
                      calc_cpu,
                      line[3]
                      ])
    print (tabulate(table, headers=headers, floatfmt='.0f'))
    print ('Total time = ' + str(tot_cpu/ 3600) + ' CPUh')
    print ()
else:
    warn('sacct table does not exist')

# If patprofile table exists, also check more detailed profiling information
if 'patprofile' in tables:
    # Only show the top x functions
    slowest_x_functs = 10
    slow_functs = []
    # First try to determine the x slowest functions
    query = db.execute('''SELECT DISTINCT Group_Function
                       FROM patprofile
                       WHERE Level>1 AND Samp_percent>5
                       ORDER BY Jobnumber ASC, Samp_percent DESC''')
    for line in islice(query, slowest_x_functs):
        slow_functs.append(line[0])
    # Next determine the x slowest lines
    slow_functs_lines = []
    query = db.execute('''SELECT DISTINCT Group_Function_Caller
                       FROM patprofile_lines
                       WHERE Level>1 AND Samp_percent>5
                       ORDER BY Jobnumber ASC, Samp_percent DESC''')
    for line in islice(query, slowest_x_functs):
        slow_functs_lines.append(line[0])
else:
    warn('patprofile table does not exist')



query = db.execute('''SELECT DISTINCT profile.Jobnumber, Dimxn, NumrawCPU, Pointspcore,
                             Numcores, Total, Nodes, Vcores_per_task, Tasks
                   FROM profile
                       JOIN stdout ON profile.Jobnumber=stdout.Jobnumber
                       JOIN jobdata ON profile.Jobnumber=jobdata.Jobnumber
                   ORDER BY Dimxn ASC, NumrawCPU ASC''')
headers = [x[0] for x in query.description]
print (tabulate(query, headers=headers, floatfmt='.0f'))
query = db.execute('''SELECT SUM(Total * Nodes * 24)
                   FROM profile
                       JOIN stdout ON profile.Jobnumber=stdout.Jobnumber
                       JOIN jobdata ON profile.Jobnumber=jobdata.Jobnumber
                   ''')
print ('Total time = ' + str(np.sum(query.fetchall())/60/60/1000) + ' CPUh')


# Now build our plots. This is inefficient, but easy. 
query = db.execute('SELECT DISTINCT NumrawCPU FROM profile ORDER BY NumrawCPU ASC')
rawcpus = np.array([x[0] for x in query])
query = db.execute('SELECT DISTINCT Dimxn FROM profile ORDER BY Dimxn ASC')
dimxns = np.array([x[0] for x in query])



# We are going to plot absolute and relative time, so create two axes
axes = [plt.figure().add_subplot(111)]
axes.append(plt.figure().add_subplot(111))
if 'patprofile' in tables:
    # Create one axis for slowest functions and one for slowest lines
    axes.append(plt.figure().add_subplot(111))
    axes.append(plt.figure().add_subplot(111))
markers = ('o', 'v', 's', 'p', '*', 'h', 'H', 'D', 'd')
# color scheme from http://colorbrewer2.org/?type=qualitative&scheme=Set1&n=9
colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3',
          '#ff7f00','#ffff33','#a65628','#f781bf','#999999']
# We have to save which marker is used for which dimxn for the legend
fake_markers = []
# For each dimxn we have in our database
minrawcpus = np.inf
for marker, dimxn in zip(markers, dimxns):
    # Reset our color cycle to plot everything in the right color
    for ax in axes:
        ax.set_prop_cycle(plt.cycler('color', colors))
    # Add a fake marker for our legend
    fake_markers.append(plt.Line2D((0,1),(0,0), color='k', marker=marker, label=dimxn, linestyle=''))
    # Initialize arrays to save plot data
    walltimes_mean = np.empty((0, len(labels)))
    walltimes_std = np.empty_like(walltimes_mean)
    if 'patprofile' in tables:
        functimes_mean = np.empty((0, len(slow_functs)))
        functimes_std = np.empty_like(functimes_mean)
        functimes_lines_mean = np.empty((0, len(slow_functs_lines)))
        functimes_lines_std = np.empty_like(functimes_lines_mean)
    # Find all different possibilities of number of raw CPUS 
    query = db.execute('''SELECT DISTINCT NumrawCPU FROM profile 
                       WHERE Dimxn=?
                       ORDER BY Pointspcore ASC''',
                       (int(dimxn),))
    
    rawcpus = np.array([x[0] for x in query])
    # For each each ppc with a specific dimxn we have in our database
    for rawcpu in rawcpus:
        # Find the walltime profiling parameters
        query = db.execute('''SELECT (Total), First_MPI_AllReduce, (Eigenmodes), (Saturation),
                               Second_MPI_AllReduce, (Initialization), (Output)
                           FROM stdout
                           JOIN profile ON stdout.Jobnumber=profile.Jobnumber
                           WHERE NumrawCPU=? AND Dimxn=?''',
                           (int(rawcpu), int(dimxn)))
        # Create an array to correctly handle duplicates
        walltimes = np.empty((0, len(labels)))
        for line in query:
            walltimes = np.vstack((walltimes, line[1:]))
        walltimes_std = np.vstack((walltimes_std, np.std(walltimes, axis=0)))
        walltimes_mean = np.vstack((walltimes_mean, np.mean(walltimes, axis=0)))


        # Opposed to the previous plot, here we have 'entries' rows instead of 4
        entries = walltimes.shape[0]
        if 'patprofile' in tables:
            # One line contains info for the 'entries' slowest functions
            func_line = np.atleast_2d(np.empty((entries, 0)))
            for func in slow_functs:
                query = db.execute('''SELECT Samp_percent
                                   FROM patprofile
                                       JOIN profile ON profile.Jobnumber=patprofile.Jobnumber
                                   WHERE NumrawCPU=? AND Dimxn=? AND Group_Function=?''',
                                   (int(rawcpu), int(dimxn), func))
                query_result = query.fetchall()
                func_entry = np.empty((entries, 1))
                # If we get a hit, put it in the 'func_entry', this will in the end be 
                # the equivalent of the provious 'line'
                if len(query_result) > 0:
                    for i, entry in enumerate(query_result):
                        func_entry[i] = entry[0]
                else:
                    # If we don't get a hit, fill the entry with NaNs
                    func_entry = np.full((entries, 1), np.nan)
                func_line = np.hstack((func_line, np.atleast_2d(func_entry)))
            functimes_std = np.vstack((functimes_std, np.std(func_line, axis=0)))
            functimes_mean = np.vstack((functimes_mean, np.mean(func_line, axis=0)))

            func_line = np.atleast_2d(np.empty((entries, 0)))
            for func in slow_functs_lines:
                query = db.execute('''SELECT Samp_percent
                                   FROM patprofile
                                       JOIN profile ON profile.Jobnumber=patprofile.Jobnumber
                                   WHERE NumrawCPU=? AND Dimxn=? AND Group_Function=?''',
                                   (int(rawcpu), int(dimxn), func))

                query_result = query.fetchall()
                func_entry = np.empty((entries, 1))
                if len(query_result) > 0:
                    for i, entry in enumerate(query_result):
                        func_entry[i] = entry[0]
                else:
                    func_entry = np.full((entries, 1), np.nan)
                func_line = np.hstack((func_line, np.atleast_2d(func_entry)))
            functimes_lines_std = np.vstack((functimes_lines_std, np.std(func_line, axis=0)))
            functimes_lines_mean = np.vstack((functimes_lines_mean, np.mean(func_line, axis=0)))

    for label, line, yerr in zip(labels, walltimes_mean.T, walltimes_std.T):
        axes[0].set_xscale('log')
        axes[0].set_yscale('log')
        axes[0].errorbar(rawcpus, line, yerr=yerr, marker=marker)
        axes[1].set_xscale('log')
        axes[1].errorbar(rawcpus, 100*line/np.sum(walltimes_mean, axis=1),
                         yerr=100*yerr/np.sum(walltimes_mean, axis=1), marker=marker)
    if 'patprofile' in tables:
        for label, line, yerr in zip(slow_functs, functimes_mean.T, functimes_std.T):
            axes[2].set_xscale('log')
            axes[2].errorbar(rawcpus, line, yerr=yerr, marker=marker)
        for label, line, yerr in zip(slow_functs_lines,
                                     functimes_lines_mean.T, functimes_lines_std.T):
            axes[3].set_xscale('log')
            axes[3].errorbar(rawcpus, line, yerr=yerr, marker=marker)

fake_colors = []
for color, label in zip(colors, labels):
    fake_colors.append(plt.Line2D((0,1),(0,0), label=label, color=color))

#Create legend from custom artist/label lists
for ax in axes[:2]:
    ax.legend(fake_colors + fake_markers,
              [marker.get_label() for marker in fake_colors] +
              ['dimxn = ' + str(marker.get_label()) for marker in fake_markers])

axes[0].set_ylabel('Walltime [ms]')
axes[0].set_xlabel('Raw CPUs')

axes[1].set_ylabel('Walltime [%]')
axes[1].set_xlabel('Raw CPUs')

if 'patprofile' in tables:
    fake_colors = []
    for color, label in zip(colors, slow_functs):
        fake_colors.append(plt.Line2D((0,1),(0,0), label=label, color=color))
    
    #Create legend from custom artist/label lists
    axes[2].legend(fake_colors + fake_markers,
                  [marker.get_label() for marker in fake_colors] +
                  ['dimxn = ' + str(marker.get_label()) for marker in fake_markers])

    fake_colors = []
    for color, label in zip(colors, slow_functs_lines):
        fake_colors.append(plt.Line2D((0,1),(0,0), label=label, color=color))
    
    #Create legend from custom artist/label lists
    axes[3].legend(fake_colors + fake_markers,
                  [marker.get_label() for marker in fake_colors] +
                  ['dimxn = ' + str(marker.get_label()) for marker in fake_markers])

    axes[2].set_ylabel('Walltime [%]')
    axes[2].set_xlabel('Raw CPUs')
    axes[3].set_ylabel('Walltime [%]')
    axes[3].set_xlabel('Raw CPUs')

plt.show()

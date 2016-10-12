import os
import subprocess
import sqlite3
import pickle
import json
from warnings import warn
import datetime

from .tabulate.tabulate import tabulate

from .qualikizrun import QuaLiKizRun, QuaLiKizBatch
from .basicpoll import database_exists


def poll_sacct(batch):
    """ Poll a job in a specific directory using sacct
    We exploit the fact that every 'run' should be unique, e.g.
    having a unique job-id.
    Arguments:
        - path: Path to poll. Path should contain a metadata
                and job file.

    Returns:
        - header: The names of the fields in the table, returned by
                  sacct
        - table:  A table containing the polled data. This table should always
                  have one row. The aforementioned headers defines the columns
                  of the table
    """
    batchdir = os.path.join(batch.batchsdir, batch.name)
    batchinfopath = os.path.join(batchdir, QuaLiKizBatch.batchinfofile)
    with open(batchinfopath, 'r') as file_:
        batchinfo = json.load(file_)

    extra_fields = 'submit,CPUTime,CPUTimeRAW,NNodes'
    cmd = 'sacct -o ' + extra_fields + ' -lXP --noconvert -j ' + batchinfo['jobnumber']
    output = subprocess.check_output(cmd, shell=True).decode('UTF-8')

    lines = output.splitlines()
    header = lines[0].split('|')
    table = [line.split('|') for line in lines[1:]]
    if len(table) > 1:
        raise Exception('Could not uniquely identify job ' + job.batch.name)
    for entry in table:
        jobname_index = header.index('JobName')
        state_index = header.index('State')
        if entry[jobname_index] == batch.name:
            if entry[state_index] != 'COMPLETED':
                warn('State of ' + entry[jobname_index] + ' is ' + entry[state_index] + ', not COMPLETED. Results might not be reliable')

    return header, table

def header_to_sql(header):
    " Helper function to generate the SQL column names """
    maxchars = max([len(x) for x in header])
    table_row = ''
    table_quest = ''
    for name in header:
        table_line = name.ljust(maxchars) + ' TEXT,'
        table_row += name + ', '
        table_quest += '?, '
        print (table_line)
    print (table_row)
    print (table_quest)

def create_sacct_database(batchlist, database_path, append=None, overwrite=None):
    """ Create a database with sacct data
    Arguments:
        -  path:          The path to be polled
        - database_path: Path to the database to be created

    Keyword Arguments:
        - overwrite: Overwrite database if exists? Default 'ask user'
        - append:    Append to table if exists? Default 'ask user'
    """
    create_table = database_exists(database_path, 'sacct', append=append, overwrite=overwrite)

    db = sqlite3.connect(database_path)
    if create_table:
        db.execute('''CREATE TABLE sacct (
                   Submit           TEXT,
                   CPUTime          TEXT,
                   CPUTimeRAW       INTEGER,
                   NNodes           INTEGER,
                   JobID            TEXT,
                   JobIDRaw         TEXT,
                   JobName          TEXT,
                   Partition        TEXT,
                   MaxVMSize        TEXT,
                   MaxVMSizeNode    TEXT,
                   MaxVMSizeTask    INTEGER,
                   AveVMSize        TEXT,
                   MaxRSS           INTEGER,
                   MaxRSSNode       TEXT,
                   MaxRSSTask       INTEGER,
                   AveRSS           INTEGER,
                   MaxPages         INTEGER,
                   MaxPagesNode     TEXT,
                   MaxPagesTask     INTEGER,
                   AvePages         INTEGER,
                   MinCPU           TEXT,
                   MinCPUNode       TEXT,
                   MinCPUTask       INTEGER,
                   AveCPU           TEXT,
                   NTasks           INTEGER,
                   AllocCPUS        INTEGER,
                   Elapsed          TEXT,
                   State            TEXT,
                   ExitCode         TEXT,
                   AveCPUFreq       TEXT,
                   ReqCPUFreqMin    TEXT,
                   ReqCPUFreqMax    TEXT,
                   ReqCPUFreqGov    TEXT,
                   ReqMem           TEXT,
                   ConsumedEnergy   INTEGER,
                   MaxDiskRead      INTEGER,
                   MaxDiskReadNode  TEXT,
                   MaxDiskReadTask  INTEGER,
                   AveDiskRead      INTEGER,
                   MaxDiskWrite     INTEGER,
                   MaxDiskWriteNode TEXT,
                   MaxDiskWriteTask INTEGER,
                   AveDiskWrite     INTEGER,
                   AllocGRES        TEXT,
                   ReqGRES          TEXT,
                   ReqTRES          TEXT,
                   AllocTRES        TEXT
                   )''')
    result = []
    for batch in batchlist:
        __, data = poll_sacct(batch) 
        db.execute('''
        INSERT INTO sacct (
        Submit, CPUTime, CPUTimeRAW, NNodes, JobID, JobIDRaw, JobName, Partition, MaxVMSize, MaxVMSizeNode, MaxVMSizeTask, AveVMSize, MaxRSS, MaxRSSNode, MaxRSSTask, AveRSS, MaxPages, MaxPagesNode, MaxPagesTask, AvePages, MinCPU, MinCPUNode, MinCPUTask, AveCPU, NTasks, AllocCPUS, Elapsed, State, ExitCode, AveCPUFreq, ReqCPUFreqMin, ReqCPUFreqMax, ReqCPUFreqGov, ReqMem, ConsumedEnergy, MaxDiskRead, MaxDiskReadNode, MaxDiskReadTask, AveDiskRead, MaxDiskWrite, MaxDiskWriteNode, MaxDiskWriteTask, AveDiskWrite, AllocGRES, ReqGRES, ReqTRES, AllocTRES)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''', data[0])
    db.commit()

    query = db.execute('select JobID, CPUTime, NNodes, AllocCPUS, Elapsed, State from sacct')
    headers = [x[0] for x in query.description]
    print (tabulate(query, headers=headers))
        
def acctstr_to_timedelta(acctstr):
    """ Covert a timestring generated with acct to timedelta """
    acctstr_split = acctstr.split(':')
    first_split = acctstr_split[0].split('-')
    minutes, seconds = [int(x) for x in acctstr_split[-2:]]
    if len(acctstr_split) < 3:
        timedelta = datetime.timedelta(seconds=seconds, minutes=minutes)
    else:
        hours = int(first_split[-1])
        if len(first_split) == 1:
            timedelta = datetime.timedelta(seconds=seconds, minutes=minutes, hours=hours)
        elif len(first_split) == 2:
            days = int(first_split[0])
            timedelta = datetime.timedelta(seconds=seconds, minutes=minutes, hours=hours, days=days)
    return timedelta

def create_database(path, database_path, append=None, overwrite=None):
    batchlist = QuaLiKizBatch.from_dir_recursive(path)
    create_sacct_database(batchlist, database_path, append=append, overwrite=overwrite)

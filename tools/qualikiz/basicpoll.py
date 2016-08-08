import os
import sqlite3
import pickle
import csv
import re

from .qualikizrun import EmptyJob, recursive_function
from .tabulate.tabulate import tabulate


def database_exists(database_path, table_name, append=None, overwrite=None):
    """ Check if database and table exists and handle collisions
    Arguments:
        - database_path: Path to the database to be checked
        - table_name:    Name of the table to be checked

    Keyword Arguments:
        - overwrite: Overwrite database if exists? Default 'ask user'
        - append:    Append to table if exists? Default 'ask user'

    Returns:
        - create_table: Does the table has to be (re-)created? Default False
    """
    if os.path.isfile(database_path):
        if overwrite is None:
            resp = input('Database exists, overwrite? [y/N]')
            if resp == '' or resp == 'N' or resp == 'n':
                overwrite = False
            elif resp == 'Y' or resp == 'y':
                overwrite = True
        if overwrite:
            print ('overwriting')
            os.remove(database_path)

    db = sqlite3.connect(database_path)
    query = db.execute("SELECT name FROM sqlite_master WHERE type='table'")
    rows = query.fetchall()
    tables = [x[0] for x in rows]

    if table_name in tables:
        if append is None:
            resp = input('Table exists, append? [Y/n]')
            if resp == '' or resp == 'Y' or resp == 'y':
                append = True
            elif resp == 'N' or resp == 'n':
                append = False
        if append:
            create_table = False
        else:
            db.execute("DROP TABLE " + table_name)
            create_table = True
    else:
        create_table = True

    return create_table

def poll_qualikiz_out(path):
    """ Read QuaLiKiz STDOUT
    This function polls for the QuaLiKiz STDOUT. It needs the metadata
    containing the jobdata. This is always generated if QuaLiKiz was
    run with one of the pythontools scripts. For example, with
    'pythontools.py inputgo runs/mini'

    Arguments:
        - path: Path to poll. Path should contain a metadata
                and job file.

    Returns:
        - header: The line containing information about the part
                  being profiled. This probably has to be parsed
                  further.
        - list:   A list with the values of the aforementioned header
    """
    with open(os.path.join(path, EmptyJob.metadatafile), 'r') as file_:
        reader = csv.reader(file_)
        job_meta = {line[0]: line[1] for line in reader}
    with open(os.path.join(path, EmptyJob.jobdatafile), 'rb') as file_:
        job = pickle.load(file_)
    profile_lines = []
    with open(os.path.join(path, job.batch.stdout), 'r') as file_:
        for line in file_:
            if line.startswith('Profiling'):
                profile_lines.append(line)

    header = ['jobnumber']
    list = [int(job_meta['jobnumber'])]
    for line in profile_lines[:3]:
        header.append(line.split('=')[0])
        num = [int(n) for n in re.split(r'\D', line) if n is not '']
        list.append(num)

    for line in profile_lines[3:]:
        header.append(line.split('=')[0])
        sec, msec = [int(n) for n in re.split(r'\D', line) if n is not '']
        list.append((sec, msec))

    # Sanity check
    if len(header) != len(list):
        raise Exception('Header and list do not have the same length')

    return header, list

def poll_jobdata(path):
    """ Read QuaLiKiz jobdata
    This is always generated if QuaLiKiz was run with one of the
    pythontools scripts. For example, with 'pythontools.py inputgo runs/mini'

    Arguments:
        - path: Path to poll. Path should contain a metadata
               and job file.
    """
    with open(os.path.join(path, EmptyJob.metadatafile), 'r') as file_:
        reader = csv.reader(file_)
        job_meta = {line[0]: line[1] for line in reader}
    with open(os.path.join(path, EmptyJob.jobdatafile), 'rb') as file_:
        job = pickle.load(file_)
    header = ['jobnumber', 'nodes', 'vcores_per_task', 'tasks']
    table = [job_meta['jobnumber'],
             job.batch.nodes,
             job.batch.vcores_per_task,
             job.batch.srun.tasks]
    table = [int(x) for x in table]
    return header, table

def create_jobdata_database(path, database_path, append=None, overwrite=None):
    """ Create a database with QuaLiKiz metadata
    Arguments:
        - path:          The path to be polled
        - database_path: Path to the database to be created

    Keyword Arguments:
        - overwrite: Overwrite database if exists? Default 'ask user'
        - append:    Append to table if exists? Default 'ask user'
    """
    create_table = database_exists(database_path, 'jobdata', append=append, overwrite=overwrite)

    db = sqlite3.connect(database_path)
    if create_table:
        db.execute('''CREATE TABLE jobdata (
            Jobnumber       INTEGER,
            Nodes           INTEGER,
            Vcores_per_task INTEGER,
            Tasks           INTEGER
          )''')
    result = recursive_function(path, poll_jobdata)

    for __, row in result:
        db.execute('''
                   INSERT INTO jobdata(
                   Jobnumber, Nodes, Vcores_per_task, Tasks)
                   VALUES (?, ?, ?, ?)''',
                   row)
    db.commit()

    query = db.execute('select * from jobdata')
    headers = [x[0] for x in query.description]
    print (tabulate(query, headers=headers, floatfmt='.0f'))

def create_qualikiz_out_database(path, database_path, append=None, overwrite=None):
    """ Create a database with parsed QuaLiKiz STDOUT
    Arguments:
        - path:          The path to be polled
        - database_path: Path to the database to be created

    Keyword Arguments:
        - overwrite: Overwrite database if exists? Default 'ask user'
        - append:    Append to table if exists? Default 'ask user'
    """
    create_table = database_exists(database_path, 'stdout', append=append, overwrite=overwrite)

    db = sqlite3.connect(database_path)
    if create_table:
        db.execute('''CREATE TABLE stdout (
            Jobnumber      INTEGER,
            Numcores       INTEGER,
            Dimx           INTEGER,
            Dimn           INTEGER,
            Eigenmodes     INTEGER,
            Saturation     INTEGER,
            Initialization INTEGER,
            Output         INTEGER,
            Total          INTEGER
          )''')
    result = recursive_function(path, poll_qualikiz_out)

    for __, row in result:
        data = [row[0]]
        data.extend(row[1])
        data.extend(row[2])
        data.extend(row[3])
        msec = [1000 * tup[0] + tup[1] for tup in row[4:]]
        data.extend(msec)

        db.execute('''
                   INSERT INTO stdout (
                   Jobnumber, Numcores, Dimx, Dimn, Eigenmodes, Saturation, Initialization, Output, Total)
                   VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)''',
                   data)
    db.commit()

    query = db.execute('select * from stdout')
    headers = [x[0] for x in query.description]
    print (tabulate(query, headers=headers, floatfmt='.0f'))

def create_database(path, database_path, append=None, overwrite=None):
    """ Create a database with basic QuaLiKiz data
    Arguments:
        - database_path: Path to the database to be created

    Keyword Arguments:
        - overwrite: Overwrite database if exists? Default 'ask user'
        - append:    Append to table if exists? Default 'ask user'
    """
    create_qualikiz_out_database(path, database_path, append=None, overwrite=None)
    create_jobdata_database(path, database_path, append=None, overwrite=False)

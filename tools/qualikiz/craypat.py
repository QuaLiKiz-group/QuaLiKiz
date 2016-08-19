import os
import subprocess
import sqlite3
import csv
from warnings import warn
from itertools import islice

import numpy as np
from .qualikizrun import recursive_function, EmptyJob
from .basicpoll import database_exists
from .tabulate.tabulate import tabulate

def profile_dir(path):
    """ Read QuaLiKiz STDOUT
    This function polls for the QuaLiKiz profiling information
    generated with CrayPat. It needs the metadata
    containing the jobdata. This is always generated if QuaLiKiz was
    run with one of the pythontools scripts. For example, with
    'pythontools.py inputgo runs/mini'

    Arguments:
        - path: Path to poll. Path should contain a metadata
                and job file.

    Returns:
        - A pair containing (jobnumber, tables), with:
            - jobnumber: the jobnumber as read from the metadata
            - All tables read from craypat. Right now the first table
              contains the sampling of line-numbers and the second table
              contains the sampling of functions
    """
    pat_result = None
    for file_ in os.listdir(path):
        if file_.startswith('QuaLiKiz+pat+'):
            pat_result = file_
    if pat_result is None:
        warn('No folder starting with QuaLiKiz+pat+ found in ' + path + ', skipping..')
    else:
        with open(os.path.join(path, EmptyJob.metadatafile), 'r') as file_:
            reader = csv.reader(file_)
            job_meta = {line[0]: line[1] for line in reader}

        cmd = 'pat_report -O ca+src,profile -s show_data=csv ' + os.path.join(path, pat_result)
        print (cmd)
        output = subprocess.check_output(cmd, shell=True, stderr=None).decode('UTF-8')
        # With the previous command we created two tables and some useless data.
        # We need to only read only the relevant parts
        read = [False, False]
        tables = [[], []]
        for line in output.splitlines():
            # Stop reading if you find an empty row
            if read[0]:
                if line == '':
                    read[0] = False
                else:
                    tables[0].append(line)
            if read[1]:
                if line == '':
                    read[1] = False
                else:
                    tables[1].append(line)
            # Start reading if you encounter the table header
            if line == 'Level,Samp%,Samp,Group/Function/Caller/PE=HIDE/Thread=HIDE':
                read[0] = True
            if line == 'Level,Samp%,Samp,Imb..Samp,Imb..Samp%,Group/Function/PE=HIDE/Thread=HIDE':
                read[1] = True

        formatted_tables = [[], []]
        # Strip the % signs from the entries
        for i, table in enumerate(tables):
            for row in table:
                formatted_tables[i].append([x.rstrip('%') for x in row.split(',')])
        result = (job_meta['jobnumber'], formatted_tables)
        return result

def create_database(path, database_path, append=None, overwrite=None):
    """ Create a database with QuaLiKiz metadata
    Arguments:
        - path:          The path to be polled
        - database_path: Path to the database to be created

    Keyword Arguments:
        - overwrite: Overwrite database if exists? Default 'ask user'
        - append:    Append to table if exists? Default 'ask user'
    """
    create_table = database_exists(database_path, 'patprofile', append=append, overwrite=overwrite)

    db = sqlite3.connect(database_path)
    if create_table:
        db.execute('''CREATE TABLE patprofile (
                Jobnumber INTEGER,
                Level INTEGER,
                Samp_percent REAL,
                Samp INTEGER,
                Imb_Samp INTEGER,
                Imb_Samp_percent REAL,
                Group_Function TEXT)''')
        db.execute('''CREATE TABLE patprofile_lines (
                Jobnumber INTEGER,
                Level INTEGER,
                Samp_percent REAL,
                Samp INTEGER,
                Group_Function_Caller TEXT)''')
    results = recursive_function(path, profile_dir)
    for result in results:
        jobnumber = result[0]
        tables = result[1]
        for i, table in enumerate(tables):
            for row in table:
                data = [jobnumber]
                data.extend(row)
                if i == 0:
                    db.execute('''
                    INSERT INTO patprofile_lines (Jobnumber, Level, Samp_percent, Samp, Group_Function_Caller)
                            VALUES (?,?,?,?,?)''', data)
                if i == 1:
                    db.execute('''
                    INSERT INTO patprofile (Jobnumber, Level, Samp_percent, Samp, Imb_Samp, Imb_Samp_percent, Group_Function)
                            VALUES (?,?,?,?,?,?,?)''', data)
    db.commit()

    query = db.execute('SELECT Level, Samp_percent, Group_Function FROM patprofile WHERE Level>1 ORDER BY Samp_percent DESC')
    headers = [x[0] for x in query.description]
    result = query.fetchmany(10)
    print (tabulate(result, headers=headers, floatfmt='.0f'))

    query = db.execute('SELECT Level, Samp_percent, Group_Function_Caller FROM patprofile_lines WHERE Level>1 ORDER BY Samp_percent DESC')
    headers = [x[0] for x in query.description]
    result = query.fetchmany(10)
    print (tabulate(result, headers=headers, floatfmt='.0f'))

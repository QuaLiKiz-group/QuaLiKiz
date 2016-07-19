import os
import datetime
import pickle
import subprocess
import sys
import math
import warnings
from warnings import warn
from collections import OrderedDict
import shutil

import numpy as np

from .edisonbatch import Srun, Batch
from . import inputfiles
from .inputfiles import QuaLiKizRun, Ion, IonList, Electron

warnings.simplefilter('always', UserWarning)
class PathException(Exception):
    def __init__(self, path):
        message = path + ' must be an absolute path or bad things will happen!'
        super().__init__(message)

class Run:
    """ A collection of QuaLiKiz Jobs
    This class is used to define a collection of QuaLiKiz Jobs, or in
    other words a collection of sbatch scripts and their input.

    Class variables:
        - scriptname:  The default name of the sbatch scripts
    """
    # pylint: disable=too-few-public-methods

    scriptname = 'edison.sbatch'
    def __init__(self, rootdir, qualikizdir="", runsdir=""):
        """ Initialize a run
        Arguments:
            - rootdir: Directory of the QuaLiKiz root repo.

        Keyword arguments:
            - runsdir:     Directory of the QuaLiKiz runs dir. Runs will be
                           created in this folder. Usually derived from
                           the rootdir
            - qualikizdir: Directory of the QuaLiKiz python tools. Usually
                           derived from the rootdir
        """
        self.jobs = {}
        if os.path.isabs(rootdir):
            self.rootdir = rootdir
        else:
            raise PathException('rootdir')

        if qualikizdir == "":
            self.qualikizdir = os.path.join(rootdir, 'tools/qualikiz')
        else:
            self.qualikizdir = qualikizdir

        if runsdir == "":
            self.runsdir = os.path.join(rootdir, 'runs')
        else:
            self.runsdir = runsdir

        if not os.path.isabs(self.qualikizdir):
            raise PathException('qualikizdir')

        if not os.path.isabs(self.runsdir):
            raise PathException('runsdir')

    def to_file(self, overwrite=False):
        """ Writes all Run folders to file. """
        # pylint: disable=invalid-name, unused-variable
        for __, job in self.jobs.items():
            job.to_file(self.runsdir, self.qualikizdir, overwrite=overwrite)


class EmptyJob:
    jobdatafile = 'jobdata.pkl'
    """ Defines an empty QuaLiKiz Job, so no input files are generated"""
    # pylint: disable=too-few-public-methods, too-many-arguments

    def __init__(self, rootdir, name, binary_basename, batch,
                 parameters_script=""):
        """ Initialize an empty QuaLiKiz job
        Arguments:
            - rootdir:         Path of the QuaLiKiz root. Usually get this
                               from Run instance
            - name:            The name of the QuaLiKiz Job. This will be the
                               name of the folder that will be generated
            - binary_basename: The name of the binary that needs to be run.
                               Relative to the batch script file
            - batch:           An instance of the Batch script

        Keyword arguments:
            - parameters_script: A string with the script used to 
                                 create the input files

        Undefined at initialization time:
            - jobnumber:  The jobnumber given by the SLURM system
            - submittime: The time the job was submitted to the SLURM system
        """
        if os.path.isabs(rootdir):
            self.binary_path = os.path.abspath(os.path.join(rootdir,
                                                        binary_basename))
        else:
            raise PathException('rootdir')

        self.parameters_script = parameters_script
        self.name = name
        self.batch = batch
        self.jobnumber = None
        self.submittime = None

    def to_file(self, runsdir, qualikizdir, overwrite=False, metapickle=True):
        """ Write all Job folders to file
        This generates a folder for each Job. Each Job has its own sbatch script.
        This will also generate all folders and links needed for QuaLiKiz to run.

        Arguments:
            - runsdir:     Directory to which to write the Jobs
            - qualikizdir: Directory of the qualikiz python module

        Keyword arguments:
            - overwrite:  Overwrite the directory if it exists. Default=False
            - metapickle: Dump the Job data to file. This is needed if you
                          want to recreate or analyze the job later. Default=True
        """
        if not os.path.isabs(qualikizdir):
            raise PathException('qualikizdir')

        if not os.path.isabs(runsdir):
            raise PathException('runsdir')

        if not os.path.isdir(runsdir):
            os.mkdir(runsdir)

        path = os.path.abspath(os.path.join(runsdir, self.name))
        try:
            os.mkdir(path)
        except FileExistsError:
            if not overwrite:
                resp = input('folder exists, overwrite? [Y/n]')
                if resp == '' or resp == 'Y' or resp == 'y':
                    overwrite = True
            if overwrite:
                print ('overwriting')
                shutil.rmtree(path)
                os.mkdir(path)
            else:
                print ('aborting')
                return 1
        if not os.path.exists(self.binary_path):
            warn('Warning! Binary at ' + self.binary_path + ' does not ' +
                 'exist! Run will fail!')
        os.symlink(self.binary_path,
                   os.path.join(path, os.path.basename(self.binary_path)))
        os.symlink(qualikizdir,
                   os.path.join(path, 'qualikiz'))
        os.makedirs(os.path.join(path, 'output/primitive'), exist_ok=True)
        os.makedirs(os.path.join(path, 'debug'), exist_ok=True)
        self.batch.to_file(path=os.path.join(path, Run.scriptname))
        if self.parameters_script == "":
            with open(os.path.join(qualikizdir,
                                   'parameters_template.py')) as file_:
                parameters_script = file_.readlines()
        else:
            parameters_script = self.parameters_script
        with open(os.path.join(path, 'parameters.py'), 'w') as file_:
              file_.writelines(parameters_script)

        if metapickle:
            with open(os.path.join(path, 'jobdata.pkl'), 'wb') as file_:
                pickle.dump(self, file_)

def generate_inputs(path, dotprint=False):
    """ Recursively create input for all jobs in a specific directory """
    if Run.scriptname in os.listdir(path):
        generate_input(path)
    else:
        for folder in os.listdir(path):
            job_folder = os.path.join(path, folder)
            if os.path.isdir(job_folder):
                if Run.scriptname in os.listdir(job_folder):
                    generate_input(job_folder)
                    print ('.', end='', flush=True)

def generate_input(job_folder):
    """ Generate input for a job in a specific directory """
    cmd = ['python', os.path.abspath(os.path.join(job_folder, 'parameters.py'))]
    subprocess.check_call(cmd)

def run_jobs(path):
    """ Recursively run all jobs in a specific directory """
    if Run.scriptname in os.listdir(path):
        run_job(path)
    else:
        for folder in os.listdir(path):
            job_folder = os.path.join(path, folder)
            if os.path.isdir(job_folder):
                if Run.scriptname in os.listdir(job_folder):
                    run_job(job_folder)

def run_job(job_folder):
    """ Run a job in a specific directory """
    input_binary = os.path.join(job_folder, 'input/p1.bin')
    if not os.path.exists(input_binary):
        warn('Warning! Input binary at ' + input_binary + ' does not ' +
             'exist! Run will fail! Please generate input binaries with ' +
             'parameters.py')


    cmd = 'git describe --tags'
    qualikiz_version = subprocess.check_output(cmd, shell=True)

    cmd = 'sbatch --workdir=' + job_folder + ' ' + os.path.join(job_folder, Run.scriptname)
    output = subprocess.check_output(cmd, shell=True)
    jobnumber = output.split()[-1]
    
    with open(os.path.join(job_folder, EmptyJob.jobdatafile), 'rb') as file_:
        job = pickle.load(file_)
        job.jobnumber = jobnumber
        job.submittime = datetime.datetime.now()
        job.qualikiz_version = qualikiz_version
    with open(os.path.join(job_folder, EmptyJob.jobdatafile), 'wb') as file_:
        pickle.dump(job, file_)
        print (jobnumber)

def poll_dirs(path):
    """ Recursively poll all jobs in a specific directory using sacct
    Arguments:
        - path: Path to recursively poll for sacct data

    Returns:
        - List of dicts containing sacct data
    """
    acct_data = []
    if EmptyJob.jobdatafile in os.listdir(path):
        acct_data.append(poll_dir(path))
    else:
        for folder in os.listdir(path):
            job_folder = os.path.join(path, folder)
            if os.path.isdir(job_folder):
                if EmptyJob.jobdatafile in os.listdir(job_folder):
                    acct_data.append(poll_dir(job_folder))
    return acct_data

def poll_dir(path):
    """ Poll a job in a specific directory using sacct
    Arguments:
        - path: Path to poll for sacct data

    Returns:
        - Dict containing sacct data
    """
    with open(os.path.join(path, EmptyJob.jobdatafile), 'rb') as file_:
        job = pickle.load(file_)
        acct_data = poll_job(job)
    return acct_data

def poll_job(job):
    """ Poll a job using sacct
    Arguments:
        - job: An EmptyJob instance. Job should have ran and have metadata

    Returns:
        - Dict containing sacct data
    """
    extra_fields = 'submit,CPUTime,CPUTimeRAW,NNodes'
    if job.jobnumber is None:
        raise Exception('Job has not run yet, cannot analyse')
    cmd = 'sacct -o ' + extra_fields + ' -lP --name ' + job.batch.name + ' -j ' + job.jobnumber.decode("utf-8")
    output = subprocess.check_output(cmd, shell=True).decode('UTF-8')
    lines = output.splitlines()
    table = []
    header = lines[0]
    headers = header.split(sep='|')
    for line in lines[1:]:
        words = line.split(sep='|')
        table.append({header: word for header, word in zip(headers, words)})
    if len(table) > 3:
        raise Exception('Could not uniquely identify job ' + job.batch.name)
    for entry in table:
        if entry['JobName'] == job.batch.name:
            if entry['State'] != 'COMPLETED':
                warn('State of ' + entry['JobName'] + ' is ' + entry['State'] + ', not COMPLETED. Results might not be reliable')
            elapsed = acctstr_to_timedelta(entry['Elapsed'])


            raw_cores = job.batch.cores_per_node * int(entry['NNodes'])
            raw_cpus = elapsed.seconds * raw_cores
    return table

def poll_to_file(acct_data, path):
    """ Writes data from poll_job(s) to file
    Arguments:
        acct_data: Data returned from poll_job(s)
        path: file to write to
    """
    if not path.endswith('.pkl'):
        warn('It is advised to add the \'.pkl\' extention,' +
             'as the file will be pickled')
    with open(path, 'wb') as file_:
        pickle.dump(acct_data, file_)


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

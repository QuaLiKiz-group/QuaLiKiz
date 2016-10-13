"""
Copyright Dutch Institute for Fundamental Energy Research (2016)
Contributors: Karel van de Plassche (karelvandeplassche@gmail.com)
License: CeCILL v2.1
"""
from warnings import warn
from qualikiz.qualikizrun import vcores_per_task
cores_per_node = 24


class Sbatch:
    """ Defines a batch job
    This class uses the OpenMP/MPI parameters as defined by Edison,
    but could in principle be extented to support more machines.

    class variables:
        - attr:             All possible attributes as defined by Edison
        - sbatch:           Names of attributes as they are in the sbatch file
        - shell:            The shell to use for sbatch scripts. Usually bash
        - repo:             The default repo to bill hours to. Usually empty
    """
    # pylint: disable=too-many-instance-attributes
    attr = ['nodes',
            'maxtime',
            'partition',
            'tasks_per_node',
            'vcores_per_task',
            'filesystem',
            'name',
            'repo',
            'stderr',
            'stdout',
            'mailtype',
            'mailuser',
            'qos']
    sbatch = ['nodes',
              'time',
              'partition',
              'ntasks-per-node',
              'cpus-per-task',
              'license',
              'job-name',
              'account',
              'error',
              'output',
              'mail-type',
              'mail-user',
              'qos']
    shell = '/bin/bash'
    repo = None
    mailtype = None
    mailuser = None
    workdir = None
    default_stderr = 'stderr.batch'
    default_stdout = 'stdout.batch'

    def __init__(self, srun_instances, name, tasks, maxtime,
                 stdout=default_stdout, stderr=default_stderr,
                 filesystem='SCRATCH', partition='regular',
                 qos='normal', HT=True):
        """ Initialize Edison batch job
        Arguments:
            - srun_instances: List of Srun instances included in the Sbatch job
            - name:           Name of the Sbatch job
            - tasks:          Amount of MPI tasks
            - maxtime:        Maximum walltime needed

        Keyword Arguments:
            - stdout:     File to write stdout to. By default 'stdout.batch'
            - stderr:     File to write stderr to. By default 'stderr.batch'
            - filesystem: The default filesystem to use. Usually SCRATCH
            - partition:  Partition to run on, for example 'debug'. By default
                          'regular'
            - qos:        Priority in the queue. By default 'normal'
            - HT:         Hyperthreading on/off. Default=True


        Calculated:
            - threads_per_core: amount of OMP threads per physical core
            - threads_per_node: amount of OMP threads per compute node
            - sockets_per_node: Amount of sockets in one compute node
            - cores_per_socket: Amount of physical CPU cores in one socket
            - cores_per_node:   Amount of physical CPU cores in one node
        """

        self.filesystem = filesystem
        if HT:
            vcores_per_core = 2  # Per definition
        else:
            vcores_per_core = 1
        # HT 48 or no HT 24
        self.vcores_per_node = cores_per_node * vcores_per_core
        self.vcores_per_task = vcores_per_task  # 2, as per QuaLiKiz
        self.tasks_per_node = int(self.vcores_per_node / self.vcores_per_task)

        nodes, remainder = divmod(tasks, self.tasks_per_node)
        self.nodes = int(nodes)
        if remainder != 0:
            self.nodes += 1
            warn(str(tasks) + ' tasks not evenly divisible over ' +
                 str(self.tasks_per_node) + ' tasks per node. Using ' +
                 str(self.nodes) + ' nodes.')

        self.qos = qos
        self.maxtime = maxtime
        self.partition = partition
        self.name = name
        self.srun_instances = srun_instances
        self.stdout = stdout
        self.stderr = stderr

    def to_file(self, path):
        """ Writes sbatch script to file
        Arguments:
            - path: Path of the sbatch script file.
        """
        sbatch_lines = ['#!' + self.shell + ' -l\n']
        for attr, sbatch in zip(self.attr, self.sbatch):
            value = getattr(self, attr)
            if value is not None:
                line = '#SBATCH --' + sbatch + '=' + str(value) + '\n'
                sbatch_lines.append(line)

        sbatch_lines.append('\nexport OMP_NUM_THREADS=' +
                            str(self.vcores_per_task) + '\n')

        # Write sruns to file
        sbatch_lines.append(self.srun_instances[0].to_string())
        for run_instance in self.srun_instances[1:]:
            sbatch_lines.append(' &\n' + run_instance.to_string())
        sbatch_lines.append('\n')

        with open(path, 'w') as file:
            file.writelines(sbatch_lines)

    @classmethod
    def from_file(cls, path):
        """ Reconstruct sbatch from sbatch file """
        new = Sbatch.__new__(cls)
        srun_strings = []
        with open(path, 'r') as file:
            for line in file:
                if line.startswith('#SBATCH --'):
                    line = line.lstrip('#SBATCH --')
                    name, value = line.split('=')
                    value = str_to_number(value.strip())
                    if name in cls.sbatch:
                        setattr(new, cls.attr[cls.sbatch.index(name)], value)
                if line.startswith('srun'):
                    srun_strings.append(line)

        new.vcores_per_node = new.tasks_per_node * new.vcores_per_task

        new.srun_instances = []
        for srun_string in srun_strings:
            new.srun_instances.append(Srun.from_string(srun_string))
        return new

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        return NotImplemented

    def __ne__(self, other):
        if isinstance(other, self.__class__):
            return not self == other
        return NotImplemented


def str_to_number(string):
    """ Convert a string in a float or int if possible """
    try:
        value = float(string)
    except ValueError:
        value = string
    else:
        if value.is_integer:
            value = int(value)
    return value


class Srun:
    """ Defines the srun command """
    default_stderr = 'stderr.run'
    default_stdout = 'stdout.run'

    def __init__(self, binary_name, tasks,
                 chdir='.',
                 stdout=default_stdout, stderr=default_stderr):
        """ Initializes the Srun class
        Arguments:
            - binary_name: The name of the binary relative to where
                           the sbatch script will be
            - tasks: Amount of MPI tasks needed for the job
        Keyword Arguments:
            - chdir:  Dir to change to before running the command
            - stdout: Standard target of redirect of STDOUT
            - stderr: Standard taget of redirect of STDERR
        """
        self.binary_name = binary_name
        self.tasks = tasks
        self.chdir = chdir
        self.stdout = stdout
        self.stderr = stderr

    def to_string(self):
        """ Create the srun string """
        string = 'srun'
        string += ' -n ' + str(self.tasks)
        string += ' --chdir ' + self.chdir
        string += ' --output ' + self.stdout
        string += ' --error ' + self.stderr
        string += ' ' + self.binary_name
        return string

    @classmethod
    def from_string(cls, string):
        """ Reconstruct the Srun from a string """
        new = Srun.__new__(cls)
        split = string.split(' ')
        new.tasks = int(split[2])
        new.chdir = split[4]
        new.stdout = split[6]
        new.stderr = split[8]
        new.binary_name = split[9].strip()
        return new

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        return NotImplemented

    def __ne__(self, other):
        if isinstance(other, self.__class__):
            return not self == other
        return NotImplemented

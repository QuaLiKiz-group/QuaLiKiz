from warnings import warn
class Batch:
    """ Defines a batch job
    This class uses the OpenMP/MPI parameters as defined by Edison,
    but could in principle be extented to support more machines.

    class variables:
        - sockets_per_node: Amount of sockets in one compute node
        - cores_per_socket: Amount of physical CPU cores in one socket
        - cores_per_node:   Amount of physical CPU cores in one compute node
        - attr:             All possible attributes as defined by Edison
        - sbatch:           Names of attributes as they are in the sbatch file
        - shell:            The shell to use for sbatch scripts. Usually bash
        - filesystem:       The default filesystem to use. Usually SCRATCH
        - repo:             The default repo to bill hours to. Usually empty
        - stdout:           The file to write stdout to. Default=qualikiz.out
        - stderr:           The file to write stderr to. Default=qualikiz.err

    functions:
        - optimize_sbatch: Calculate tasks_per_node and nodes based on HT
        - optimize_sbatch_general: LEGACY
        - to_file: Write edison sbatch to file
    """
    # pylint: disable=too-many-instance-attributes
    sockets_per_node = 2
    cores_per_socket = 12
    cores_per_node = sockets_per_node * cores_per_socket
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
    filesystem = 'SCRATCH'
    repo = None
    stdout = 'qualikiz.out'
    stderr = 'qualikiz.err'
    mailtype = None
    mailuser = None
    workdir = None

    def __init__(self, partition, srun, qos='normal', HT=True):
        """ Initialize Edison batch job
        Arguments:
            - partition: Partition to run on, for example 'debug'
            - srun:      Instance of the Srun class; describes run parameters

        Keyword arguments:
            - qos:       Priority in the queue. Default='normal'
            - HT:        Hyperthreading on/off. Default=True

        Undefined at initialization time:
            - name:             Name of the sbatch job
            - nodes:            Amount of nodes to run on
            - maxime:           Maximum allowed walltime
            - tasks_per_node:   Amount of MPI tasks per node
            - threads_per_task: Amount of OMP threads per MPI task
            - vcores_per_task:   Amount of virtual cores per MPI task

        Calculated:
            - threads_per_core: amount of OMP threads per physical core
            - threads_per_node: amount of OMP threads per compute node
        """
        self.qos = qos
        self.nodes = None #-N --nodes
        self.maxtime = None
        self.partition = partition
        self.tasks_per_node = None #-n, --ntasks-per-node
        self.threads_per_task = None #OMP_NUM_THREADS=2
        self.name = None
        self.srun = srun
        self.vcores_per_task = None #-c --cpus-per-task

        if HT:
            self.threads_per_core = 2
        else:
            self.threads_per_core = 1
        self.threads_per_node = (self.sockets_per_node *
                                 self.cores_per_socket *
                                 self.threads_per_core)

    def optimize_sbatch(self):
        """ Calculate optimal sbatch parameters
        Some of the parameter validation is kept in place for legacy purposes,
        even though they should be always constant as long as QuaLiKiz
        does not change.

        Constants:
            - self.threads_per_task = 2, as per QuaLiKiz code
            - self.vcores_per_task   = 2, max 1 MPI task per virtual core

        Calculated:
            - self.tasks_per_node: Amount of MPI tasks per compute node. Should be 12
                                   without HT and 24 with HT
            - self.nodes:          Amount of nodes used for batch job. It should be at least
                                   one node.
        """
        self.threads_per_task = 2 # Stuck as per QuaLiKiz code
        threads_per_vcore = 1 # Never give one (virtual) CPU more than one task
        self.vcores_per_task, remainder = divmod(self.threads_per_task, threads_per_vcore)
        if remainder != 0:
            warn(str(self.threads_per_task) + ' threads per task not evenly divisible over ' + \
                 str(self.threads_per_core) + ' threads per cpu. Using ' + \
                 str(self.vcores_per_task) + ' cpus per tasks.')

        # 12 without HT, 24 with HT
        self.tasks_per_node, remainder = divmod(self.threads_per_node, self.threads_per_task)
        if remainder != 0:
            warn(str(self.threads_per_node) + ' threads per node not evenly divisible over ' + \
                 str(self.threads_per_task) + ' threads per task. Using ' + \
                 str(self.tasks_per_node) + ' tasks per node.')
        nodes, remainder = divmod(self.srun.tasks, self.tasks_per_node)
        if nodes == 0:
            nodes = 1
            warn('Cannot use 0 nodes for ' + str(self.srun.tasks) + ' tasks. Using 1 node instead.')
        else:
            if remainder != 0:
                nodes += 1
                warn(str(self.srun.tasks) + ' tasks not evenly divisible over ' + \
                     str(self.tasks_per_node) + ' tasks per node. Using ' + \
                     str(nodes) + ' nodes.')

        self.nodes = nodes
        if self.threads_per_task * self.tasks_per_node > self.threads_per_node:
            raise Exception('Too many cpus per node!')

    def optimize_sbatch_general(self, tasks_per_node):
        """ LEGACY """
        raise Exception('Function kept in for legacy purpose! Use at your own risk!')
        # pylint: disable=unreachable
        if tasks_per_node > 48:
            raise Exception('Too many tasks per node!')
        self.tasks_per_node = tasks_per_node # Usually 1-4, max 24 or 48 with HT
        self.threads_per_task, remainder = divmod(self.threads_per_node, self.tasks_per_node)
        if remainder != 0:
            self.threads_per_task += 1
            warn(str(self.tasks_per_node) + ' tasks per node not evenly divisible over ' + \
                 str(self.threads_per_node) + ' threads per node. Using ' + \
                 str(self.threads_per_task) + ' threads per task.')
        threads_per_cpu = 1
        self.vcores_per_task, remainder = divmod(self.threads_per_task, threads_per_cpu)
        if remainder != 0:
            warn(str(self.threads_per_task) + ' threads per task not evenly divisible over ' + \
                 str(self.threads_per_core) + ' threads per cpu. Using ' + \
                 str(self.vcores_per_task) + ' cpus per tasks.')

        self.nodes, remainder = divmod(self.srun.tasks, tasks_per_node)
        if remainder != 0:
            self.nodes += 1
            warn(str(self.srun.tasks) + ' tasks not evenly divisible over ' + \
                 str(self.tasks_per_node) + ' tasks per node. Using ' + \
                 str(self.nodes) + ' nodes.')

        if self.threads_per_task * self.tasks_per_node > self.threads_per_node:
            raise Exception('Too many cpus per node!')

    def to_file(self, path='edison.sbatch'):
        """ Writes sbatch script to file
        Remember to set all attributes defined in the __init__ function, or the batch
        script will not work

        Keyword arguments:
            - path: (relative) path of the sbatch script file. This will be the name
                    of the file. Default=edison.sbatch
        """
        sbatch_lines = ['#!' + self.shell + ' -l\n']
        for attr, sbatch in zip(self.attr, self.sbatch):
            value = getattr(self, attr)
            if  value != None:
                line = '#SBATCH --' + sbatch + '=' + str(value) + '\n'
                sbatch_lines.append(line)

        sbatch_lines.append('\nexport OMP_NUM_THREADS=' + str(self.threads_per_task) + '\n')
        sbatch_lines.append('srun -n ' + str(self.srun.tasks) + ' ' +self.srun.binary_name + '\n')

        with open(path, 'w') as file:
            file.writelines(sbatch_lines)


class Srun:
    """ Defines the srun command inside a sbatch script"""
    # pylint: disable=too-few-public-methods
    def __init__(self, binary_name, tasks):
        """ Initializes the Srun class
        Arguments:
            - binary_name: The name of the binary relative to where the sbatch script will be
            - tasks: Amount of MPI tasks needed for the job
        """
        self.binary_name = binary_name
        self.tasks = tasks

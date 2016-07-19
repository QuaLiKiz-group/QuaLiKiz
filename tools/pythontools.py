""" Command line interface for Python tools """
import sys
import os
from qualikiz import compare
from qualikiz import qualikizrun
from qualikiz import legacy
import subprocess
import os, inspect


if len(sys.argv) < 2:
    raise Exception('Please supply a command')
command = sys.argv[1]

if command == 'compare':
    if len(sys.argv) < 4:
        raise Exception('Please supply two paths to compare')
    path1 = sys.argv[2]
    path2 = sys.argv[3]
    compare.compare_runs(path1, path2)

elif command == 'dump':
    if len(sys.argv) < 4:
        raise Exception('Please supply the style and path of the file')
    style = sys.argv[2]
    path = sys.argv[3]
    if style == 'bin':
        array = compare.bin_to_np(path)
    elif style == 'ascii':
        array = compare.ascii_to_np(path)
    else:
        raise Exception('Unknown file style \'' + style + '\'')
    print (array)

elif command == 'poll':
    if len(sys.argv) < 3:
        raise Exception('Please supply poll path')
    if len(sys.argv) == 4:
        targetdir = sys.argv[3]
    else:
        targetdir = './polldump.pkl'
    path = sys.argv[2]
    acctdata = qualikizrun.recursive_function(path, qualikizrun.poll_dir)
    qualikizrun.poll_to_file(acctdata, targetdir)

elif command == 'inputgo':
    if len(sys.argv) < 3:
        raise Exception('Please supply run path')
    path = sys.argv[2]
    qualikizrun.recursive_function(path, qualikizrun.generate_input)
    qualikizrun.run_jobs(path)

elif command == 'convertto':
    if len(sys.argv) < 3:
        raise Exception('Please supply current style')
    current = sys.argv[2]
    if len(sys.argv) < 4:
        raise Exception('Please supply conversion target')
    target = sys.argv[3]
    if len(sys.argv) < 5:
        raise Exception('Please supply input path to convert')
    path = sys.argv[4]
    if current == 'current':
        if target == '2.3.0':
            legacy.convert_current_to_2_3_0(path)
        else:
            raise Exception('Unkown target')
    else:
        raise Exception('Unkown current style')

else:
    testsdir =  os.path.abspath(
        os.path.join(os.path.abspath(inspect.getfile(inspect.currentframe())),
                     '../../tests'))
    if command == 'create':
        if len(sys.argv) < 3:
            raise Exception('Please supply create target')
        create_target = sys.argv[2]
        if len(sys.argv) == 4:
            target_dir = sys.argv[3]
        else:
            target_dir = ''

        if create_target == 'mini':
            cmd = ['python',
                   os.path.join(testsdir, 'mini.py'),
                   target_dir ]
        elif create_target == 'performance':
            cmd = ['python',
                   os.path.join(testsdir, 'performance.py'),
                   target_dir]
        else:
            raise Exception('Unknown create target \'' + create_target + '\'')
        subprocess.check_call(cmd)

    elif command == 'analyse':
        if len(sys.argv) < 3:
             raise Exception('Please supply path to analyse')
        path = sys.argv[2]
        cmd = ['python',
               os.path.join(testsdir, 'performance_analyse.py'),
               path]
        subprocess.check_call(cmd)

    else:
        raise Exception('Unknown command \'' + command + '\'')

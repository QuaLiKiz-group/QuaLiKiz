#!/usr/bin/env python3
""" Command line interface for Python tools
Copyright Dutch Institute for Fundamental Energy Research (2016)
Contributors: Karel van de Plassche (karelvandeplassche@gmail.com)
License: CeCILL v2.1
"""
import sys
import os
import subprocess
import inspect

from qualikiz import (compare, legacy, qualikizrun, basicpoll,
                      sacct, craypat, fs_manipulation)


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
    print(array)

elif command == 'poll':
    if len(sys.argv) < 3:
        raise Exception('Please supply poll path')
    if len(sys.argv) == 4:
        targetdir = sys.argv[3]
    else:
        targetdir = './polldb.sqlite3'
    path = sys.argv[2]
    basicpoll.create_database(path, targetdir)
    sacct.create_database(path, targetdir, overwrite=False)
    craypat.create_database(path, targetdir, overwrite=False)

elif command == 'basicpoll':
    if len(sys.argv) < 3:
        raise Exception('Please supply poll path')
    if len(sys.argv) == 4:
        targetdir = sys.argv[3]
    else:
        targetdir = './polldb.sqlite3'
    path = sys.argv[2]
    basicpoll.create_database(path, targetdir)

elif command == 'sacctpoll':
    if len(sys.argv) < 3:
        raise Exception('Please supply poll path')
    if len(sys.argv) == 4:
        targetdir = sys.argv[3]
    else:
        targetdir = './polldb.sqlite3'
    path = sys.argv[2]
    sacct.create_database(path, targetdir)

elif command == 'patpoll':
    if len(sys.argv) < 3:
        raise Exception('Please supply poll path')
    if len(sys.argv) == 4:
        targetdir = sys.argv[3]
    else:
        targetdir = './polldb.sqlite3'
    path = sys.argv[2]
    craypat.create_database(path, targetdir)

elif command == 'movecomplete':
    if len(sys.argv) < 4:
        raise Exception('please supply source and target')
    else:
        sourcedir = sys.argv[2]
        targetdir = sys.argv[3]
    fs_manipulation.move_completed(sourcedir, targetdir)

elif command == 'inputgo':
    if len(sys.argv) < 3:
        raise Exception('Please supply run path')
    path = sys.argv[2]
    batchlist = qualikizrun.QuaLiKizBatch.from_dir_recursive(path)
    for batch in batchlist:
        batch.prepare()
        batch.generate_input()
        batch.queue_batch()

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
        if target == '2.3.1' or target == '2.3.0':
            legacy.convert_current_to_2_3_2(path)
            legacy.convert_2_3_2_to_2_3_1(path)
        elif target == '2.3.2':
            legacy.convert_current_to_2_3_2(path)
        else:
            raise Exception('Unkown target')
    elif current == '2.3.2':
        if target == '2.3.1' or target == '2.3.0':
            legacy.convert_2_3_2_to_2_3_1(path)
        else:
            raise Exception('Unkown target')
    else:
        raise Exception('Unkown current style')

else:
    testsdir = os.path.abspath(
        os.path.join(os.path.abspath(inspect.getfile(inspect.currentframe())),
                     '../tests'))
    if command == 'create':
        if len(sys.argv) < 3:
            raise Exception('Please supply create target')
        create_target = sys.argv[2]
        if len(sys.argv) == 4:
            target_dir = sys.argv[3]
        else:
            target_dir = ''

        if create_target == 'mini':
            cmd = ['python3',
                   os.path.join(testsdir, 'mini.py'),
                   target_dir]
        elif create_target == 'performance':
            cmd = ['python3',
                   os.path.join(testsdir, 'performance.py'),
                   target_dir]
        else:
            raise Exception('Unknown create target \'' + create_target + '\'')
        subprocess.check_call(cmd)

    elif command == 'analyse':
        if len(sys.argv) < 3:
            raise Exception('Please supply path to analyse')
        path = sys.argv[2]
        cmd = ['python3',
               os.path.join(testsdir, 'performance_analyse.py'),
               path]
        subprocess.check_call(cmd)

    else:
        raise Exception('Unknown command \'' + command + '\'')

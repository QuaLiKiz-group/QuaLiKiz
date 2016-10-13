#!/usr/bin/env python3
import sys
import os
sys.path.append(os.path.dirname(sys.path[0]))
from qualikiz.qualikizrun import QuaLiKizBatch, QuaLiKizRun
from qualikiz.edisonbatch import Srun, Sbatch

command = 'go'
if len(sys.argv) == 2:
    if sys.argv[1] != '':
        command = sys.argv[1]

batch = QuaLiKizBatch.from_dir('.')

if command.startswith('input'):
    batch.generate_input()

if command.endswith('go'):
    batch.queue_batch()

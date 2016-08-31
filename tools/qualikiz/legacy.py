import os
import shutil
from warnings import warn
def convert_current_to_2_3_2(inputdir):
    input_table = {'p1' : 'dimx',
                   'p2' : 'dimn',
                   'p3' : 'nions',
                   'p4' : 'phys_meth',
                   'p5' : 'coll_flag',
                   'p6' : 'rot_flag',
                   'p7' : 'verbose',
                   'p8' : 'separateflux',
                   'p9' : 'numsols',
                   'p10': 'relacc1',
                   'p11': 'relacc2',
                   'p12': 'maxruns',
                   'p13': 'maxpts',
                   'p14': 'timeout',
                   'p15': 'R0',
                   'p16': 'kthetarhos',
                   'p17': 'x',
                   'p18': 'rho',
                   'p19': 'Ro',
                   'p20': 'Rmin',
                   'p21': 'Bo',
                   'p22': 'qx',
                   'p23': 'smag',
                   'p24': 'alphax',
                   'p25': 'Machtor',
                   'p26': 'Autor',
                   'p27': 'Machpar',
                   'p28': 'Aupar',
                   'p29': 'gammaE',
                   'p30': 'Te',
                   'p31': 'ne',
                   'p32': 'Ate',
                   'p33': 'Ane',
                   'p34': 'typee',
                   'p35': 'anise',
                   'p36': 'danisdre',
                   'p37': 'Ti',
                   'p38': 'normni',
                   'p39': 'Ati',
                   'p40': 'Ani',
                   'p41': 'typei',
                   'p42': 'anisi',
                   'p43': 'danisdri',
                   'p44': 'Ai',
                   'p45': 'Zi'}
    for old_name, new_name in input_table.items():
        try:
            os.rename(os.path.join(inputdir, new_name + '.bin'),
                      os.path.join(inputdir, old_name + '.bin'))
        except FileNotFoundError:
            warn('File ' + os.path.join(inputdir, new_name + '.bin') + ' not found, skipping..')

def convert_2_3_2_to_2_3_1(inputdir):
    input_table = {'p8' : 'p9',
                   'p9' : 'p10',
                   'p10': 'p11',
                   'p11': 'p12',
                   'p12': 'p13',
                   'p13': 'p14',
                   'p14': 'p16',
                   'p15': 'p17',
                   'p16': 'p18',
                   'p17': 'p19',
                   'p18': 'p20',
                   'p19': 'p21',
                   'p20': 'p15',
                   'p21': 'p22',
                   'p22': 'p23',
                   'p23': 'p24',
                   'p24': 'p25',
                   'p25': 'p26',
                   'p26': 'p27',
                   'p27': 'p28',
                   'p28': 'p29',
                   'p29': 'p30',
                   'p30': 'p31',
                   'p31': 'p32',
                   'p32': 'p33',
                   'p33': 'p34',
                   'p34': 'p35',
                   'p35': 'p36',
                   'p36': 'p44',
                   'p37': 'p45',
                   'p38': 'p37',
                   'p39': 'p38',
                   'p40': 'p39',
                   'p41': 'p40',
                   'p42': 'p41',
                   'p43': 'p42',
                   'p44': 'p43'}

    tempdir = os.path.join(inputdir, 'temp')
    os.mkdir(tempdir)
    suffix = '.bin'
    for file in os.listdir(inputdir):
        if file.endswith(suffix):
            os.rename(os.path.join(inputdir, file), os.path.join(tempdir, file))

    for old_name, new_name in input_table.items():
        try:
            os.rename(os.path.join(tempdir, new_name + suffix),
                      os.path.join(inputdir, old_name + suffix))
        except FileNotFoundError:
            warn('File ' + os.path.join(inputdir, new_name + suffix) + ' not found, skipping..')

    os.remove(os.path.join(tempdir, 'p8.bin'))

    for file in os.listdir(tempdir):
        os.rename(os.path.join(tempdir, file), os.path.join(inputdir, file))

    os.rmdir(tempdir)

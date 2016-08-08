import numpy as np
from qualikiz.inputfiles import *
""" Imitates the old MATLAB script """
scan_names = ['iAt', 'iT', 'qx', 'smag']
xpoints = 3
scan_ranges = [np.linspace(2, 12, xpoints),
               np.linspace(2.4, 24, xpoints),
               np.linspace(1, 5, xpoints),
               np.linspace(0.1, 3, xpoints)]

elec = Electron(T=8., n=5., At=9., An=3., type=1, anis=1., danisdr=0.)
D = Ion(name='main_D', Ai=2., Zi=1., n=0.8, T=8., At=0., An=3., type=1, anis=1., danisdr=0.)
Be = Ion(name='Be', Ai=9., Zi=4., n=0.1, T=8., At=0., An=2.9, type=1, anis=1., danisdr=0.)
W = Ion(name='W+42', Ai=184., Zi=42., n=0.0, T=8., At=0., An=3., type=3, anis=1., danisdr=0.)

ions = IonList(D, Be, W)

npoints = 8
kthetarhos = np.append(np.linspace(0.1, 0.8, npoints),
                       np.linspace(6, 48, npoints))
dict_ = {
        'R_0':     3,
        'x':        .5,
        'rho':      .5,
        'Ro':      3.,
        'Rmin':    1.,
        'Bo':      3.,
        'R0':      3.,
        'qx':      2.,
        'smag':    1.,
        'alphax':  0.,
        'Machtor': 0.,
        'Autor':   0.,
        'ninorm1': True,
        'Ani1':    True,
        'QN_grad': True
        }

run = QuaLiKizRun(scan_names, scan_ranges, kthetarhos, elec, ions, **dict_)
hyperedge = run.setup_hyperedge()
run.to_file(hyperedge)


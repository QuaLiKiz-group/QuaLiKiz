#!/bin/python3
import array
import os
import sys
from os import path
from collections import OrderedDict
import itertools
import numpy as np

def to_bytelist_proto(self):
    """ Protype of converting a class to a list of bytes
    """
    bytelist = [array.array('d', [value]) for value in self.values()]
    return bytelist


class Particle(OrderedDict):
    """ Particle (ion or electron)
    """
    keynames = ['T', 'n', 'At', 'An', 'type', 'anis', 'danisdr']
    def __init__(self, **kwargs):
        """ Initialize a Particle.
        Usually it is better to create an Electron or Ion instead.

    kwargs:
        T :       Temperature in keV
        n :       Density in 10^19 m^-3 for electrons, relative factor to
                  electron denisity for ions
        At:       Normalized logarithmic temperature gradient
                  A_t = - (R/T) * (dT/dr)
        An:       Normalized logarithmic density gradient
                  A_n = - (R/n) * (dn/dr)
        type: 1:  active
              2:  adiabatic
              3:  passing at ion scales
        anise:    Temperature anisotropy T_perp / T_para at LFS
        danisedr: Radial gradient of temperature anisotropy

        kwargs (ion only):
            Ai: Ion mass in amu
            Zi: Ion charge in e
        """
        values = [kwargs[arg] for arg in Particle.keynames]
        super().__init__(zip(Particle.keynames, values))

    def to_nparray(self):
        nparray = np.fromiter(self.values(), float)
        return nparray


class Electron(Particle):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def to_bytelist(self):
        return to_bytelist_proto(self)


class Ion(Particle):
    keynames = ['Ai', 'Zi']
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.name = kwargs['name']
        self['Ai'] = kwargs['Ai']
        self['Zi'] = kwargs['Zi']


class IonList(list):
    """ Convenient wrapper for a list of Ions

        args:
            list of Ions

    """
    def __init__(self, *args):
        super().__init__(args)

    def to_bytelist(self):
        ionarray = [array.array('d', range(len(self))) for i in range(9)]
        for i, ion in enumerate(self):
            for j, entry in enumerate(ion.values()):
                ionarray[j][i] = entry

        return ionarray

    def __setitem__(self, key, value):
        for ion in self:
            ion[key] = value
        return self


class QuaLiKizRun(OrderedDict):
    """ Wraps multi-dimensional QuaLiKiz scan input files

    This class contains multiple dicts sorted by meaning. Those dicts are:
    elec: an Electron that describes the electrons in the plasma
    ions: an IonList with all ions contained in the plasma
    meta: a Meta instance with all values that don't change for different
          radial points
    special: a Special instance with all values that need special treatment
             when writing to binary
    spatial: a Spatial instance with all values that change for different
             radial points

    """
    def __init__(self, _scan_names, _scan_ranges, kthetarhos, electrons, ions, **kwargs):
        """ Initialize a single QuaLiKiz run
        Normally you would initialize your QuaLiKiz run with this, then set
        up a N-dimensional scan with setup_hypercube and write the
        QuaLiKiz input files with to_file

        args:
            _scan_names:  Names of the parameters to be scanned. Ion parameters
                         can be scanned by prepending 'i' and electron
                         parameters can be scanned by prepending 'e'. If you
                         want to scan over a specific ions parameter, you can
                         prepend 'in' for example, 'i1'.
            _scan_values: The values to be scanned matching the names in
                         _scan_names
            kthetarhos:  The wave spectrum to be scanned
            electrons:   An Electron instance describing the electron
                         population
            ions:        An IonList instance describing the ion population

        kwargs:
            all kwargs described in the Meta, Special and Spatial classes
            ninorm1:     Flag to normalize main ion concentration to maintain
                         quasineutrality
            Ani1:        Flag to normalize main ion gradient to maintain
                         quasineutrality
            QN_grad:     Flag for maintaining quasineutrality of gradients

        """
        super().__init__()
        self._scan_ranges = _scan_ranges
        self._scan_names = _scan_names

        numscan = np.product([len(range) for range in _scan_ranges])
        numwave = len(kthetarhos)
        nion = len(ions)
        self['meta'] = self.Meta(numscan=numscan, numwave=numwave, nion=nion, **kwargs)
        self['special'] = self.Special(kthetarhos)
        self['spatial'] = self.Spatial(**kwargs)
        self['elec'] = electrons
        self['ions'] = ions

        self.ninorm1 = kwargs.get('ninorm1', True)
        self.Ani1 = kwargs.get('Ani1', True)
        self.QN_grad = kwargs.get('QN_grad', True)
        self.normalize()

    def newscan(self, scan_names, scan_ranges):
        """ Set up a new scan based on the old QuaLiKizRun """
        kwargs = {}
        kwargs.update(self['meta'])
        kwargs.pop('numwave')
        kwargs.pop('numscan')
        kwargs.pop('nion')
        kwargs.update(self['special'])
        kwargs.pop('kthetarhos')
        kwargs.update(self['spatial'])
        kwargs.update({'ninorm1': self.ninorm1,
                       'Ani1': self.Ani1,
                       'QN_grad': self.QN_grad})
        return QuaLiKizRun(scan_names, scan_ranges,
                           self['special']['kthetarhos'],
                           self['elec'], self['ions'], **kwargs)
    def normalize(self):
        if self.ninorm1:
            self.normalize_density()
        if self.Ani1:
            self.normalize_gradient()
        if self.QN_grad:
            self.check_quasi()

    def normalize_density(self):
        ninormd = [ion['n'] if ion['type'] != 3 and ion['type'] != 4 else 0 for ion in self['ions']]
        if self.ninorm1 and len(self['ions']) > 1: #sets ninorm of 1st species to maintian quasineutrality
            self['ions'][0]['n'] = (1 - np.sum([ninorm * ion['Zi'] for ninorm, ion in zip(ninormd[1:], self['ions'][1:])]))/self['ions'][0]['Zi']

    def normalize_gradient(self):
        ninormd = [ion['n'] if ion['type'] != 3 and ion['type'] != 4 else 0 for ion in self['ions']]
        if self.Ani1 and len(self['ions']) > 1: #sets Ani of 1st species to maintian quasineutrality
            self['ions'][0]['An'] = (self['elec']['An'] - np.sum([ninorm * ion['An'] * ion['Zi'] for ninorm, ion in zip(ninormd[1:], self['ions'][1:])]))/(self['ions'][0]['Zi'] * self['ions'][0]['n'])

    def check_quasi(self):
        quasicheck = [ion['Zi'] * ion['n'] for ion in self['ions']]
        quasicheck_grad = [ion['Zi'] * ion['n'] * ion['An'] for ion in self['ions']]
        quasitol = 1e-5


    class Meta(OrderedDict):
        """ Wraps variables that stay constant during the whole QuaLiKiz run """
        keynames = ['numscan', 'numwave', 'nion', 'phys_meth', 'coll_flag',
                    'rot_flag', 'verbose', 'numsols', 'relacc1', 'relacc2',
                    'maxruns', 'maxpts', 'timeout', 'R_0']
        def __init__(self, **kwargs):
            """ Initialize Meta class
            kwargs:
                numscan:   Number of scan points
                numwave:   Number of wavenumbers
                nions:     Number of ions in the system
                phys_meth: Flag for additional calculation
                coll_flag: Flag for collisionality
                rot_flag:  Flag for rotation
                verbose:   Flag for level of output verbosity
                numsols:   Number of requested solutions
                relacc1:   Relative accuracy of 1D integrals
                relacc2:   Relative accuracy of 2D integrals
                maxruns:   Number of runs jumping directly to Newton between
                           contour checks
                maxpts:    Number of integrant evaluations done in 2D integral
                timeout:   Upper time limit [s] for wavenumber/scan point
                           solution finding
                R_0:       [m] Geometric major radius used for normalizations
            """
            defaults = {
                'numscan':   None,
                'numwave':   None,
                'nion':      None,
                'phys_meth': True,
                'coll_flag': True,
                'rot_flag':  False,
                'verbose':   True,
                'numsols':   3,
                'relacc1':   1e-3,
                'relacc2':   2e-2,
                'maxruns':   1,
                'maxpts':    5e5,
                'timeout':   60,
                'R_0':       None
            }
            values = [kwargs.get(arg, defaults[arg]) for arg in self.keynames]
            super().__init__(zip(self.keynames, values))

        def to_bytelist(self):
            return to_bytelist_proto(self)


    class Special(OrderedDict):
        """ Wraps variables that need special convertion to binary"""
        def __init__(self, kthetarhos):
            """ Initialize Special class
            kwargs:
                kthetarhos: Wave spectrum input
            """
            super().__init__(kthetarhos=kthetarhos)

        def to_bytelist(self):
            return [array.array('d', value) for value in self.values()]


    class Spatial(OrderedDict):
        """ Wraps variables that change per scan point """
        in_args = ['x', 'rho', 'Ro', 'Rmin', 'Bo', 'qx', 'smag',
                   'alphax', 'Machtor', 'Autor']
        extra_args = ['Machpar', 'Aupar']

        def __init__(self, **kwargs):
            """ Initialize Spatial class
            kwargs:
                x:       [m] Radial normalized coordinate
                rho:     [-] Normalized toroidal flux coordinate
                Ro:      [m] Major radius
                Rmin:    [m] Minor radius
                Bo:      [T] Magnetic field at magnetic axis
                qx:      Local q-profile value
                smag:    Local magnetic shear s def rq'/q
                alphax:  Local MHD alpha
                Machtor: Normalized toroidal velocity
                Autor:   Toroidal velocity gradient

            internally calculated:
                Machpar: Normalized parallel velocity
                Aupar:   Parallel velocity gradient
                gammaE:  Normalized perpendicular ExB flow shear
            """
            super().__init__()
            for arg in self.in_args:
                self[arg] = kwargs[arg]
            Machtor = self['Machtor']
            Autor = self['Autor']
            qx = self['qx']
            epsilon = self['x']/self['Ro']
            self['Machpar'] = Machtor / np.sqrt(1+(epsilon/qx)**2)
            self['Aupar'] = Autor/ np.sqrt(1+(epsilon/qx)**2)
            self['gammaE'] = -epsilon/qx*Autor # - Machtor./qx.*(1-smag./qx)   #ExB shear velocity normalized by cref/R. This form assumes pure toroidal rotation in the problem, but free to choose any number

        def to_bytelist(self):
            return to_bytelist_proto(self)

    def howto_isequal(self, name):
        """ Creates function that checks if a variable has a specific value
        This function hides the internal complexity of the QuaLiKizRun
        class. You can check a specific internal variable, or check for
        an Electron variable by prepending 'e', or check for all Ions
        by prepending 'i'. You can also check for a specific Ion with
        'i#', for example 'i1'
        """
        if name.startswith('i'):
            try:
                ionnumber = int(name[1])
                name = name[2:]
                def isequal(self, value):
                    return self['ions'][ionnumber][name] == value
            except ValueError:
                name = name[1:]
                def isequal(self, value):
                    return all([ion[name] == value for ion in self['ions']])

            if (not name in Ion.keynames) and (not name in Particle.keynames):
                raise Exception(name + ' does not exist!')

        else:
            if name in self.Spatial.in_args + self.Spatial.extra_args:
                def isequal(self, value):
                    return self['spatial'].__getitem__(name) == value
            elif name.startswith('e'):
                name = name[1:]
                def isequal(self, value):
                    return self['elec'].__getitem__(name) == value
        return isequal

    def howto_setitem(self, name):
        """ Created a function that can set a value by keyname
        Use this method to set a value in the QuaLiKizRun class. This is prefered
        to setting a value by hand as you would do in a normal dict.
        This is because of the multi-layered structure of the QuaLiKizRun class.
        You can set a specific internal variable, or set
        an Electron variable by prepending 'e', or set all Ions
        by prepending 'i'. You can also set a specific Ion with
        'i#', for example 'i1'
        """
        if name.startswith('i'):
            try:
                ionnumber = int(name[1])
                name = name[2:]
                if (name == 'n' or name == 'Zi') and self.ninorm1:
                    def set(self, value):
                        self['ions'][ionnumber][name] = value
                        self.normalize_density()
                elif name == 'An':
                    def set(self, value):
                        self['ions'][ionnumber][name] = value
                        self.normalize_gradient()
                else:
                    def set(self, value):
                        self['ions'][ionnumber][name] = value
            except ValueError:
                name = name[1:]
                if (name == 'n' or name == 'Zi') and self.ninorm1:
                    def set(self, value):
                        self['ions'].__setitem__(name, value)
                        self.normalize_density()
                elif name == 'An':
                    def set(self, value):
                        self['ions'].__setitem__(name, value)
                        self.normalize_gradient()
                else:
                    def set(self, value):
                        self['ions'].__setitem__(name, value)

            if (not name in Ion.keynames) and (not name in Particle.keynames):
                raise Exception(name + ' does not exist!')

        else:
            if name in self.Spatial.in_args + self.Spatial.extra_args:
                def set(self, value):
                    self['spatial'].__setitem__(name, value)
            elif name.startswith('e'):
                name = name[1:]
                if (name == 'n' or name == 'Zi') and self.ninorm1:
                    def set(self, value):
                        self['elec'].__setitem__(name, value)
                        self.normalize_density()
                elif name == 'An':
                    def set(self, value):
                        self['elec'].__setitem__(name, value)
                        self.normalize_gradient()
                else:
                    def set(self, value):
                        self['elec'].__setitem__(name, value)
        return set

    def setup_hypercube(self):
        set_items = [self.howto_setitem(name) for name in self._scan_names]
        isequal_items = [self.howto_isequal(name) for name in self._scan_names]

        numbytes = self['meta']['numscan']
        spat_bytelist = [array.array('d', [0] * numbytes) for i in range(13)]
        elec_bytelist = [array.array('d', [0] * numbytes) for i in range(7)]
        numbytes = self['meta']['numscan'] * self['meta']['nion']
        ions_bytelist = [array.array('d', [0] * numbytes) for i in range(9)]

        hypercube_scan_values = itertools.product(*self._scan_ranges)
        for numscan, scan_values in enumerate(hypercube_scan_values):
            for i, scan_value in enumerate(scan_values):
                if not isequal_items[i](self, scan_value):
                    set_items[i](self, scan_value)

            # Here we exploit the ordered nature of the spatial and elec dicts.
            # We iterate over both of them and put them in our initialized
            # array.
            values = itertools.chain(self['spatial'].values(),
                                     self['elec'].values())
            for i, value in enumerate(values):
                (spat_bytelist + elec_bytelist)[i][numscan] = value

            # Note that the ion array is in C ordering, not F ordering
            for j, ion in enumerate(self['ions']):
                for k, value in enumerate(ion.values()):
                    ions_bytelist[k][j * self['meta']['numscan'] + numscan] = value

        # Some magic because electron type is a QuaLiKizRun constant
        elec_bytelist[4] = array.array('d', [elec_bytelist[4][0]])
        return (self['meta'].to_bytelist() + self['special'].to_bytelist() +
                spat_bytelist + elec_bytelist + ions_bytelist)


    def to_file(self, bytelist):
        """ Writes a bytelist to file """
        inputdir = path.join(path.dirname(sys.argv[0]), 'input')
        try:
            os.mkdir(inputdir)
        except FileExistsError:
            pass

        for i, byte in enumerate(bytelist, 1):
            with open(path.join(inputdir, 'p' + str(i) + '.bin'), 'wb') as file:
                byte.tofile(file)

        # Magic to fix ordering of the p-files
        os.rename(path.join(inputdir, 'p14.bin'), 'temp')
        for i in range(15, 21):
            os.rename(path.join(inputdir, 'p' + str(i) + '.bin'),
                      path.join(inputdir, 'p' + str(i - 1) + '.bin'))
        os.rename('temp', path.join(inputdir, 'p20.bin'))


        os.rename(path.join(inputdir, 'p43.bin'), 'temp')
        os.rename(path.join(inputdir, 'p44.bin'), 'temp2')
        for i in reversed(range(36, 43)):
            os.rename(path.join(inputdir, 'p' + str(i) + '.bin'),
                      path.join(inputdir, 'p' + str(i + 2) + '.bin'))
        os.rename('temp', path.join(inputdir, 'p36.bin'))
        os.rename('temp2', path.join(inputdir, 'p37.bin'))

    def legacy(self):
        #
        ##Set the number and range of wave number points
        #
        #
        ##NOTE: for general scans, any of the below can be changed to a vector of size (ones(scann,1))
        ##e.g. for a q-profile scan: qx=linspace(1,4,scann)'
        #
        qe = 1.602176565e-19
        mp = 1.672621777e-27
        cref = np.sqrt(qe*1e3/mp)
        #cthi = np.sqrt(2*qe*Tix(:,1)*1e3/Ai(1)/mp)
        cthi = np.sqrt(2*qe*D['T']*1e3/D['Ai']/mp)
        #
        #
        ##Calculate auxilliary quantities for plotting
        Epsilonx=Rmin*x/Ro
        ft=2*(2*Epsilonx)**(0.5)/np.pi; #trapped particle fraction
        tau=elec['T']/D['T']
#for i=1:scann
        ninormz= [ion['n'] if ion['type'] != 3 else 0 for ion in ions]
        Ziz2 = [ion['Zi'] ** 2 for ion in ions]
        Zeffx = np.sum([x * y for x, y in zip(ninormz, Ziz2)])

def reference_example():
    """ Imitates the old MATLAB script """
    scan_names = ['iAt']
    scan_ranges = [np.linspace(6,12,3)]

    elec = Electron(T=8., n=5., At=9., An=3., type=1, anis=1., danisdr=0.)
    D = Ion(name='main_D', Ai=2., Zi=1., n=0.8, T=8., At=6., An=3., type=1, anis=1., danisdr=0.)
    Be = Ion(name='Be', Ai=9., Zi=4., n=0.1, T=8., At=6., An=2.9, type=1, anis=1., danisdr=0.)
    W = Ion(name='W+42', Ai=184., Zi=42., n=0.0, T=8., At=6., An=3., type=3, anis=1., danisdr=0.)

    ions = IonList(D, Be, W)

    kthetarhos = np.linspace(0.1, 0.8, 8)
    dict = {
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

    sim = QuaLiKizRun(scan_names, scan_ranges, kthetarhos, elec, ions, **dict)
    return sim
#   Zeffx(i)=sum(ninormz(i,:).*Zi(i,:).^2); #calculate Zeff
#end
#
#figure;
#disp('Figures displayed such that the input data can be verified by eye')
#subplot(221)
#value(gca,'FontSize',18)
#plot(1:scann,Ati(:,1),'r',1:scann,Tix(:,1),'r--',1:scann,Ate,'b',1:scann,Tex,'b--','LineWidth',2)
#xlabel('Scan index')
#
#l1=legend('-R\nabla{T_i}/T_i','Ti','-R\nabla{T_e}/T_e','T_e');
#value(l1,'FontSize',12)
#legend boxoff
#grid on
#subplot(222)
#value(gca,'FontSize',18)
#plot(1:scann,Ane,'g',1:scann,Nex,'g--','LineWidth',2)
#xlabel('Scan index')
#
#l2=legend('-R\nabla{n_e}/n_e','n_e');
#value(l2,'FontSize',12)
#legend boxoff
#grid on
#subplot(223)
#value(gca,'FontSize',18)
#plot(1:scann,tau,'c',1:scann,Zeffx,'c--',1:scann,ft,'c-.','LineWidth',2)
#xlabel('Scan index')
#
#l3=legend('T_e/T_i','Z_{eff}','f_t');
#value(l3,'FontSize',12)
#legend boxoff
#grid on
#subplot(224)
#value(gca,'FontSize',18)
#plot(1:scann,smag,'m',1:scann,qx,'m--',1:scann,alphax,'m-.s','LineWidth',2)
#xlabel('Scan index')
#l4=legend('s','q','\alpha');
#value(l4,'FontSize',12)
#legend boxoff
#grid on
#
#### End
#
#if(1)
#
#### Saving binary files
### ---------------------------------
### MANAGING QUALIKIZ INPUT VARIABLES
### ---------------------------------
if __name__ == "__main__":
    sim = reference_example()
    hypercube = sim.setup_hypercube()
    sim.to_file(hypercube)
    #
    #scan_names = ['T', 'n', 'At', 'An', 'anis', 'danisdr', 'Ai', 'Zi']
    #scan_names = ['i' + i for i in scan_names]
    #scan_ranges = [np.linspace(6,12,4) for i in scan_names]
    #
    #sim.newscan(scan_names, scan_ranges)
    #scan_names = ['iAt', 'eAt', 'eAn', 'i1n', 'iT', 'q']
    #
    filenames = sorted(os.listdir('../runs/run/input'))
    print (filenames)
    for filename in filenames:
        try:
            with open('input/' + filename, 'rb') as file:
                print (filename)
                arr = array.array('d')
                try:
                    arr.fromfile(file, 100)
                except EOFError:
                    pass
            print (arr)
        except FileNotFoundError:
            pass
        with open('../runs/run/input/' + filename, 'rb') as file:
            arr2 = array.array('d')
            try:
                arr2.fromfile(file, 100)
            except EOFError:
                pass
            print (arr2)

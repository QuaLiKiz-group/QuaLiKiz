"""
Copyright Dutch Institute for Fundamental Energy Research (2016)
Contributors: Karel van de Plassche (karelvandeplassche@gmail.com)
License: CeCILL v2.1
"""
import array
import copy
import json
import itertools
from collections import OrderedDict
from warnings import warn

import numpy as np
import scipy as sc
# import scipy.optimize


class Particle(dict):
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
            danisdr: Radial gradient of temperature anisotropy

            kwargs (ion only):
                Ai: Ion mass in amu
                Zi: Ion charge in e
        """
        values = [kwargs[arg] for arg in Particle.keynames]
        super().__init__(zip(Particle.keynames, values))


class Electron(Particle):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)


class Ion(Particle):
    keynames = ['A', 'Z']

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self['A'] = kwargs['A']
        self['Z'] = kwargs['Z']


class IonList(list):
    """ Convenient wrapper for a list of Ions

        args:
            list of Ions

    """
    def __init__(self, *args):
        super().__init__(args)

    def __setitem__(self, key, value):
        for ion in self:
            ion[key] = value
        return self


class QuaLiKizXpoint(dict):
    """ A single xpoint in a QuaLiKiz run.

    Typically a QuaLiKiz run scans over multiple xpoints. This class
    contains multiple dicts sorted by meaning. Those dicts are:
        elec: an Electron that describes the electrons in the plasma
        ions: an IonList with all ions contained in the plasma
        meta: a Meta instance with all values that don't change for different
              radial points
        special: a Special instance with all values that need special treatment
                 when writing to binary
        geometric: a Geometry instance with all values that change for
                   different radial points

    """
    def __init__(self, kthetarhos, electrons, ions, **kwargs):
        """ Initialize a single QuaLiKiz run
        Normally you would initialize your QuaLiKiz run with this, then set
        up a N-dimensional scan with setup_hypercube and write the
        QuaLiKiz input files with to_file

        args:
            kthetarhos:  The wave spectrum to be scanned
            electrons:   An Electron instance describing the electron
                         population
            ions:        An IonList instance describing the ion population

        kwargs:
            all kwargs described in the Meta, Special and Geometry classes
            ninorm1:     Flag to normalize main ion concentration to maintain
                         quasineutrality
            Ani1:        Flag to normalize main ion gradient to maintain
                         quasineutrality
            QN_grad:     Flag for maintaining quasineutrality of gradients
            x_rho:       Flag to keep rho and x equal if set with __setitem__
        """
        super().__init__()

        dimn = len(kthetarhos)
        nions = len(ions)
        self['meta'] = self.Meta(dimn=dimn, nions=nions, **kwargs)
        self['special'] = self.Special(kthetarhos)
        self['geometry'] = self.Geometry(**kwargs)
        self['elec'] = electrons
        self['ions'] = ions

        self['norm'] = {}
        self['norm']['ninorm1'] = kwargs.get('ninorm1', True)
        if self['norm']['ninorm1']:
            self.normalize_density()
        self['norm']['Ani1'] = kwargs.get('Ani1', True)
        if self['norm']['Ani1']:
            self.normalize_gradient()
        self['norm']['QN_grad'] = kwargs.get('QN_grad', True)
        if self['norm']['QN_grad']:
            self.check_quasi()
        self['norm']['x_rho'] = kwargs.get('x_rho', True)

    def normalize_density(self):
        """ Set density of 1st ion to maintian quasineutrality """
        if len(self['ions']) > 1:
            ions = filter(lambda x: x['type'] < 3, self['ions'][1:])
            self['ions'][0]['n'] = (1 -
                                    sum(ion['n'] * ion['Z'] for ion in ions) /
                                    self['ions'][0]['Z'])

    def normalize_gradient(self):
        """ Set density gradient of 1st ion to maintian quasineutrality """
        if len(self['ions']) > 1:
            ions = filter(lambda x: x['type'] < 3, self['ions'][1:])
            self['ions'][0]['An'] = ((self['elec']['An'] -
                                      sum(ion['n'] * ion['An'] * ion['Z']
                                          for ion in ions)) /
                                     (self['ions'][0]['Z'] *
                                      self['ions'][0]['n']))

    def check_quasi(self):
        """ Check if quasineutrality is maintained """
        ions = filter(lambda x: x['type'] < 3, self['ions'])
        quasicheck = abs(sum(ion['n'] * ion['Z'] for ion in ions) - 1)
        ions = filter(lambda x: x['type'] < 3, self['ions'])
        quasicheck_grad = abs(sum(ion['n'] * ion['An'] * ion['Z']
                                  for ion in ions) - self['elec']['An'])
        quasitol = 1e-5
        if quasicheck > quasitol:
            raise Exception('Quasineutrality violated!')
        if quasicheck_grad > quasitol:
            raise Exception('Quasineutrality gradient violated!')

    def match_zeff(self, zeff):
        """ Adjust ni1 to match the given Zeff """
        ions = filter(lambda x: x['type'] < 3, self['ions'][2:])
        if len(self['ions']) > 1:
            self['ions'][1]['n'] = ((zeff -
                                     self['ions'][0]['n'] *
                                     self['ions'][0]['Z'] ** 2 -
                                     sum(ion['n'] * ion['Z'] ** 2
                                         for ion in ions)) /
                                    self['ions'][1]['Z'] ** 2)
        # Sanity check
        # print (np.isclose(self.calc_zeff(), zeff))

    def calc_zeff(self):
        """ Calculate Zeff """
        ions = filter(lambda x: x['type'] < 3, self['ions'])
        return sum(ion['n'] * ion['Z'] ** 2 for ion in ions)

    def match_nustar(self, nustar):
        """ Set Te to match the given Nustar """
        # First set everything needed for nustar: Zeff, Ne, q, R0, Rmin, x
        # Rewrite formula for nustar to form nustar = c1 / Te^2 (c2 + ln(Te))
        zeff = self.calc_zeff()
        c1 = (6.9224e-5 * zeff * self['elec']['n'] * self['geometry']['qx'] *
              self['geometry']['Ro'] *
              (self['geometry']['Rmin'] * self['geometry']['x'] /
               self['geometry']['Ro']) ** -1.5)
        c2 = 15.2 - 0.5 * np.log(0.1 * self['elec']['n'])

        def Tex(x): return c1 / x ** 2 * (c2 + np.log(x)) - nustar

        # Initial guess
        Tex0 = np.sqrt(c1 * c2 / nustar)
        self['elec']['T'] = sc.optimize.newton(Tex, Tex0)

        # nustar_calc = c1 / Te ** 2 * (c2 + np.log(Te))
        # Sanity check
        # print(np.isclose(nustar_calc, nustar))

    def match_tite(self, tite):
        """ Set all Ions temperature to match the given Ti/Te """
        self['ions']['T'] = tite * self['elec']['T']

    class Meta(dict):
        """ Wraps variables that stay constant during the QuaLiKiz run """
        keynames = ['phys_meth', 'coll_flag',
                    'rot_flag', 'verbose', 'separateflux',
                    'numsols', 'relacc1', 'relacc2',
                    'maxruns', 'maxpts', 'timeout', 'R0']

        def __init__(self, **kwargs):
            """ Initialize Meta class
            kwargs:
                phys_meth: Flag for additional calculation of output parameters
                coll_flag: Flag for collisionality
                rot_flag:  Flag for rotation
                verbose:   Flag for level of output verbosity
                separateflux: Flag for toggling output of separate
                              ITG, TEM, ETG fluxes
                numsols:   Number of requested solutions
                relacc1:   Relative accuracy of 1D integrals
                relacc2:   Relative accuracy of 2D integrals
                maxruns:   Number of runs jumping directly to Newton between
                           contour checks
                maxpts:    Number of integrant evaluations done in 2D integral
                timeout:   Upper time limit [s] for wavenumber/scan point
                           solution finding
                R0:       [m] Geometric major radius used for normalizations
            """
            defaults = {
                'phys_meth': 2,
                'coll_flag': True,
                'rot_flag':  False,
                'verbose':   True,
                'separateflux':   False,
                'numsols':   3,
                'relacc1':   1e-3,
                'relacc2':   2e-2,
                'maxruns':   1,
                'maxpts':    5e5,
                'timeout':   60,
                'R0':       None
            }
            values = [kwargs.get(arg, defaults[arg]) for arg in self.keynames]
            super().__init__(zip(self.keynames, values))

    class Special(dict):
        """ Wraps variables that need special convertion to binary"""
        def __init__(self, kthetarhos):
            """ Initialize Special class
            kwargs:
                kthetarhos: Wave spectrum input
            """
            super().__init__(kthetarhos=kthetarhos)

    class Geometry(dict):
        """ Wraps variables that change per scan point """
        in_args = ['x', 'rho', 'Ro', 'Rmin', 'Bo', 'qx', 'smag',
                   'alphax', 'Machtor', 'Autor']
        extra_args = ['Machpar', 'Aupar', 'gammaE']

        def __init__(self, **kwargs):
            """ Initialize Geometry class
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
            self['Aupar'] = Autor / np.sqrt(1+(epsilon/qx)**2)
            self['gammaE'] = -epsilon/qx*Autor

    def __setitem__(self, key, value):
        """ Set value in nested dict
        Use this method to set a value in the QuaLiKizRun class.
        It adds some extra abstraction for the multi-layered structure
        of the QuaLiKizRun class. You can set a specific internal variable,
        or set an Electron variable by appending 'e', or set all Ions
        by appending 'i'. You can also set a specific Ion with
        'i#', for example 'i1'
        """
        if key.endswith('i') or key[-1].isdigit():
            if key[-1].isdigit():
                ionnumber = int(key[-1])
                key = key[:-2]
                if (key == 'n' or key == 'Z') and self['norm']['ninorm1']:
                    self['ions'][ionnumber][key] = value
                    self.normalize_density()
                elif key == 'An' and self['norm']['Ani1']:
                    self['ions'][ionnumber][key] = value
                    self.normalize_gradient()
                else:
                    self['ions'][ionnumber][key] = value
                if ((key == 'n' or key == 'Z' or key == 'An')
                   and self['norm']['QN_grad']):
                    self.check_quasi()
            else:
                key = key[:-1]
                if (key == 'n' or key == 'Zi') and self['norm']['ninorm1']:
                    self['ions'].__setitem__(key, value)
                    self.normalize_density()
                elif key == 'An' and self['norm']['Ani1']:
                    self['ions'].__setitem__(key, value)
                    self.normalize_gradient()
                else:
                    self['ions'].__setitem__(key, value)
                if ((key == 'n' or key == 'Z' or key == 'An')
                   and self['norm']['QN_grad']):
                    self.check_quasi()

            if (key not in Ion.keynames) and (key not in Particle.keynames):
                raise Exception(key + ' does not exist!')
        elif key.endswith('e'):
            key = key[:-1]
            if (key == 'n' or key == 'Zi') and self['norm']['ninorm1']:
                self['elec'].__setitem__(key, value)
                self.normalize_density()
            elif key == 'An' and self['norm']['Ani1']:
                self['elec'].__setitem__(key, value)
                self.normalize_gradient()
            else:
                self['elec'].__setitem__(key, value)
            if ((key == 'n' or key == 'Z' or key == 'An')
               and self['norm']['QN_grad']):
                self.check_quasi()
        elif key in self.Geometry.in_args + self.Geometry.extra_args:
            if key == 'x' and self['norm']['x_rho']:
                self['geometry'].__setitem__('rho', value)
            if key == 'rho' and self['norm']['x_rho']:
                self['geometry'].__setitem__('x', value)
            self['geometry'].__setitem__(key, value)
        elif key == 'Zeff':
            self.match_zeff(value)
        elif key == 'Nustar':
            self.match_nustar(value)
        elif key == 'Ti_Te_rel':
            self.match_tite(value)
        else:
            super().__setitem__(key, value)


class QuaLiKizPlan(dict):
    """ Wrapper for the scan plan for a QuaLiKiz run

    This class can be used to define a scan plan, in other words, over which
    values will be scanned in the QuaLiKiz run. This is given in the form of
    a xpoint base and a strategy how the exact points will be generated
    from this base. Usually this is a line or its N-D equivalent the edges of
    a hyperrectangle, or a hyperrectangle itself.
    """
    def __init__(self, scan_dict, scan_type, xpoint_base):
        """ Initialize the QuaLiKizPlan

        args:
            scan_dict:   Dictionary with as keys the names of the variables to
                         be scanned and as values the values to be scanned.
                         Use an OrderedDict to conserve ordering.
            scan_type:   How the points are generated. Currently accepts
                         'hyperedge' and 'hyperrect'
            xpoint_base: The xpoint base used as base for the generation
        """
        self['scan_dict'] = scan_dict
        self['scan_type'] = scan_type
        self['xpoint_base'] = xpoint_base

    def calculate_dimx(self):
        """ Calculate the amount of xpoints, also known as dimx

        This depends on the scan_type
        """
        lenlist = [len(range) for range in self['scan_dict'].values()]
        if self['scan_type'] == 'hyperedge':
            dimx = int(np.sum(lenlist))
        elif self['scan_type'] == 'hyperrect':
            dimx = int(np.product(lenlist))
        else:
            raise Exception('Unknown scan_type \'' + self['scan_type'] + '\'')
        return dimx

    def calculate_dimxn(self):
        """ Calculate dimxn
        """
        kthetarhos = self['xpoint_base']['special']['kthetarhos']
        return self.calculate_dimx() * len(kthetarhos)

    def edge_generator(self):
        """ Generates the points on the edge of a hyperrectangle
        """
        intersec = [x[0] for x in self['scan_dict'].values()]
        # yield intersec
        for i, (__, values) in enumerate(self['scan_dict'].items()):
            for value in values:
                point = copy.deepcopy(intersec)
                point[i] = value
                # if point != intersec:
                yield point

    def setup(self):
        """ Set up the QuaLiKiz scan

        Pass the binary generator the correct generator depending on the
        scan_type
        """
        if self['scan_type'] == 'hyperedge':
            names = list(self['scan_dict'].keys())
            bytes = self.setup_scan(names, self.edge_generator())
        elif self['scan_type'] == 'hyperrect':
            values = itertools.product(*self['scan_dict'].values())
            names = list(self['scan_dict'].keys())
            bytes = self.setup_scan(names, values)
        else:
            raise Exception('Unknown scan_type \'' + self['scan_type'] + '\'')
        return bytes

    def _sanity_check_setup(self, scan_names):
        """ Check if the order of scan_names is correct """
        try:
            index = scan_names.index('Zeff')
        except ValueError:
            pass
        else:
            if any(name in scan_names[index:] for name in ['ni', 'ni0']):
                warn('Warning! Set ni before setting Zeff')
        try:
            index = scan_names.index('Nustar')
        except ValueError:
            pass
        else:
            if any(name in scan_names[index:] for name in
                   ['Zeff', 'ne', 'q', 'Ro', 'Rmin', 'x', 'rho']):
                warn('Warning! Set Zeff, ne, q, Ro, Rmin, x' +
                     'and rho before setting Nustar')
        try:
            index = scan_names.index('Ti_Te_rel')
        except ValueError:
            pass
        else:
            if any(name in scan_names[index:] for name in ['Te']):
                warn('Warning! Set Te before setting Ti_Te_rel')

    def setup_scan(self, scan_names, scan_list):
        """ Set up a QuaLiKiz scan

        scan_names should be the names of the parameters being scanned over.
        This is a list with the same length of list-like objects generated
        by scan_list. Scan_list should be a generator (or list of lists)
        that generates the values matching the values of the scan_names.
        """
        self._sanity_check_setup(scan_names)
        dimxpoint = copy.deepcopy(self['xpoint_base'])

        # Initialize all the arrays that will eventually be written to file
        dimx = self.calculate_dimx()
        dimn = len(dimxpoint['special']['kthetarhos'])
        nions = len(dimxpoint['ions'])

        bytes = dict(zip(QuaLiKizXpoint.Geometry.in_args +
                         QuaLiKizXpoint.Geometry.extra_args,
                         [array.array('d', [0] * dimx) for i in range(13)]))
        bytes.update(dict(zip([x + 'e' for x in Electron.keynames],
                              [array.array('d', [0] * dimx)
                               for i in range(7)])))
        dimxi = dimx * nions
        bytes.update(dict(zip([x + 'i' for x in Electron.keynames +
                               Ion.keynames],
                              [array.array('d', [0] * dimxi)
                               for i in range(9)])))
        calc = {'dimx': dimx,
                'dimn': dimn,
                'nions': nions}

        # Put the three numbers we already calculated in an array
        for key, value in calc.items():
            bytes[key] = array.array('d', [value])

        # Rename what we call 'ni' to what QuaLiKiz calls 'normni'
        bytes['normni'] = bytes.pop('ni')

        numscan = -1
        # Iterate over the scan_list, each next() should provide a list-like
        # object with as many entries as we have different parameters
        for scan_values in scan_list:
            numscan += 1
            # Set the dimxn point value to the value in the list.
            for scan_name, scan_value in zip(scan_names, scan_values):
                dimxpoint[scan_name] = scan_value

            # Now iterate over all the values in the xpoint dict and add them
            # to our array
            for name, value in dimxpoint['geometry'].items():
                bytes[name][numscan] = value
            for name, value in dimxpoint['elec'].items():
                bytes[name + 'e'][numscan] = value

            # Note that the ion array is in C ordering, not F ordering
            for j, ion in enumerate(dimxpoint['ions']):
                for name, value in ion.items():
                    if name == 'n':
                        name = 'normn'
                    bytes[name + 'i'][j * dimx + numscan] = value

        # Some magic because electron type is a QuaLiKizRun constant
        bytes['typee'] = array.array('d', [bytes['typee'][0]])

        for name, value in dimxpoint['special'].items():
            bytes[name] = array.array('d', value)

        for name, value in dimxpoint['meta'].items():
            bytes[name] = array.array('d', [value])
        return bytes

    def to_json(self, filename):
        """ Dump the QuaLiKiz plan to json file

        The QuaLiKiz plan, including the xpoint base, can be fully
        recontructed later using the from_json function
        """
        with open(filename, 'w') as file_:
            json.dump(self, file_, indent=4)

    @classmethod
    def from_json(cls, filename):
        """ Load the QuaLiKiz plan from json

        Reconstruct the QuaLiKiz plan based on the given json file.
        Backwards compatibility is not guaranteed, so preferably
        generate the json with the same version as which you load it
        with.
        """
        with open(filename, 'r') as file_:
            data = json.load(file_, object_pairs_hook=OrderedDict)
            scan_dict = data.pop('scan_dict')
            scan_type = data.pop('scan_type')

            kthetarhos = data['xpoint_base']['special'].pop('kthetarhos')
            data['xpoint_base'].pop('special')
            elec = Electron(**data['xpoint_base'].pop('elec'))
            ionlist = []
            for ion in data['xpoint_base']['ions']:
                ionlist.append(Ion(**ion))
            ions = IonList(*ionlist)
            data['xpoint_base'].pop('ions')
            dict_ = {}
            for dicts in data['xpoint_base'].values():
                dict_.update(dicts)

            xpoint_base = QuaLiKizXpoint(kthetarhos, elec, ions, **dict_)
            return QuaLiKizPlan(scan_dict, scan_type, xpoint_base)

#!/bin/python3
import numpy as np
import struct
import array
import os
class Particle:
    def __init__(self, **kwargs):
        self.T = kwargs['T']
        self.At = kwargs['At']
        self.An = kwargs['An']
        self.type = kwargs['type']
        self.anis = kwargs['anis']
        self.danisdr = kwargs['danisdr']
        #prop_defaults = {"T": None, # keV
        #                 "At": None, #logarithmic gradient normalized to R
        #                 "An": None, # arbitrary for main ions, will be rewritten to ensure quasineutrality of gradients
        #                 "type": None, #1 for active, 2 for adiabatic, 3 for tracer (also won't be included in quasineutrality checks), 4 for tracer but included in Zeff
        #                 "anis": None, #Tperp/Tpar at LFS. Tix defined as Tpar
        #                 "danisdr": None # d/dr(Tperp/Tpar) at LFS. Tix defined as Tpar
        #                 }
        #for (prop, default) in prop_defaults.items():
        #    setattr(self, prop, kwargs.get(prop, default))


class Electron(Particle):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.n = kwargs['n']
        #prop_defaults = {
        #                 "n": None, #ne is in units of 10^19
        #                 }

        #for (prop, default) in prop_defaults.items():
        #    setattr(self, prop, kwargs.get(prop, default))


class Ion(Particle):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.name = kwargs['name']
        self.Ai = kwargs['Ai']
        self.Zi = kwargs['Zi']
        self.ninorm = kwargs['ninorm']
        #prop_defaults = {"name": None,
        #                 "Ai": None,
        #                 "Zi": None,
        #                 "ninorm": None, #ni/ne arbitrary for main ions (will be rewritten to ensure quasineutrality)
        #                 }

        #for (prop, default) in prop_defaults.items():
        #    setattr(self, prop, kwargs.get(prop, default))
from collections import OrderedDict
class Sim(OrderedDict):
    def __init__(self):
        scan_name = 'iAt'
        scan_value = np.linspace(6,12,3)
        self['elec'] = elec = Electron(T=8., n=5., At=9., An=3., type=1, anis=1., danisdr=0.)
        D = Ion(name='main_D', Ai=2., Zi=1., ninorm=0.8, T=8., At=6., An=3., type=1, anis=1., danisdr=0.)
        Be = Ion(name='Be', Ai=9., Zi=4., ninorm=0.1, T=8., At=6., An=2.9, type=1, anis=1., danisdr=0.)
        W = Ion(name='W+42', Ai=184., Zi=42., ninorm=0.0, T=8., At=6., An=3., type=3, anis=1., danisdr=0.)
        self['ions'] = ions = [D, Be, W]
        self['meta'] = OrderedDict()
        self['meta']['numscan'] = len(scan_value)
        kthetarhos = np.linspace(0.1,0.8,8)
        self['meta']['numwave'] = len(kthetarhos)
        self['meta']['nion'] = len(ions)
        self['meta']['phys_meth'] = 1.0# Method on additional integrals on particle flux if not 0
        self['meta']['coll_flag'] = 1.0# Flag for collisionless simulations if 0
        self['meta']['rot_flag']  = 0.0# Do not include rotation if 0
        self['meta']['verbose']   = 1  # Default verbose output
        self['meta']['numsols']   = 3  # Number of solutions in the output
        nprocs = 8 # Number of processors in parallel computation
        self['meta']['relacc1']   = 1e-3 #relative accuracy in 1D integrals
        self['meta']['relacc2']  = 2e-2 #relative accuracy in 2D integrals
        self['meta']['maxruns']   = 1 #number of runs between contour checks
        self['meta']['maxpts']    = 5e5 #number of integrand evaluations done in 2D integral
        self['meta']['timeout'] = 60 #seconds after which further solution searching is stopped, at a given wavenumber and radius
        self['meta']['kthetarhos'] = kthetarhos
        self['meta']['x']= x = 0.5	# x is the r/a throughout the scan: impacts the fraction of trapped particles	
        self['meta']['rho'] = self['meta']['x'] #assumes circular geometry
        self['meta']['Ro'] = Ro = 3. #major radius
        self['meta']['Rmin'] = Rmin = 1. #minor radius
        self['meta']['Bo'] = Bo = 3. #magnetic field
        self['meta']['R0'] = R0 = 3 #geometric major radius
        self['meta']['qx'] = qx = 2.00 #set q-profile
        self['meta']['smag'] =1.00 #set magnetic shear
        self['meta']['alphax']=0. #MHD alpha
        #rotation inputs. Intrinsic assumption that all species rotate together

        self['meta']['Machtor'] = Machtor = 0.0 # Toroidal mach number: vtor / cref . Can use cthi value above to scale to desired value with respect to main ion
        self['meta']['Autor'] = Autor = 0.0 # Toroidal velocity shear: dvtor/dr * R/cref 
        epsilon = x/Ro
        self['meta']['Machpar'] = Machtor / np.sqrt(1+(epsilon/qx)**2)   # Parallel mach number: vpar / cref . This form assumes pure toroidal rotation in the problem, but free to choose any number
        self['meta']['Aupar'] = Autor/ np.sqrt(1+(epsilon/qx)**2)  #Parallel velocity shear: vpar/dr * R/cref . This form assumes pure toroidal rotation in the problem, but free to choose any number

        self['meta']['gammaE'] = -epsilon/qx*Autor # - Machtor./qx.*(1-smag./qx)   #ExB shear velocity normalized by cref/R. This form assumes pure toroidal rotation in the problem, but free to choose any number
        ntheta=64 #resolution of parallel direction
        set_ninorm1=1      #flag for automating main ion concentration for maintaining QN
        set_Ani1=1 #flag for automating main ion gradient for maintaining QN
        set_QN_grad=1 #flag for maintaining quasineutrality of gradients.
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
        cthi = np.sqrt(2*qe*D.T*1e3/D.Ai/mp)
        #
        #
        ##Calculate auxilliary quantities for plotting
        Epsilonx=Rmin*x/Ro
        ft=2*(2*Epsilonx)**(0.5)/np.pi; #trapped particle fraction
        tau=elec.T/D.T
        #
        #Quasineutrality (including gradients) test and Zeff calculation
        ninormd= [ion.ninorm if ion.type != 3 and ion.type != 4 else 0 for ion in ions] 
        if set_ninorm1 and len(ions) > 1: #sets ninorm of 1st species to maintian quasineutrality
            ions[0].ninorm = (1 - np.sum([ninorm * ion.Zi for ninorm, ion in zip(ninormd[1:], ions[1:])]))/ions[0].Zi
        
        if set_Ani1 and len(ions) > 1: #sets Ani of 1st species to maintian quasineutrality
            ions[0].An = (elec.An - np.sum([ninorm * ion.An * ion.Zi for ninorm, ion in zip(ninormd[1:], ions[1:])]))/(ions[0].Zi * ions[0].ninorm)
        
        quasicheck = [ion.Zi * ion.ninorm for ion in ions]
        quasicheck_grad = [ion.Zi * ion.ninorm * ion.An for ion in ions]
        quasitol = 1e-5
        print (np.sum(quasicheck) - 1 < quasitol)
        print (np.sum(quasicheck_grad) - elec.An < quasitol)
        #    error(['Quasineutrality not respected by an absolute value ',num2str(abs(quasicheck(i))),'. Tolerance is ',num2str(quasitol),'. Scan index = ',num2str(i)]);
        #  error(['Quasineutrality of gradients not respected by an absolute value ',num2str(abs(quasicheck_grad(find(noQNgrad==1))')),'. Tolerance is ',num2str(quasitol),'. Scan indices = ',num2str(find(noQNgrad==1)')]);
        #end
        #
        ninormz= [ion.ninorm if ion.type != 3 else 0 for ion in ions] 
        Ziz2 = [ion.Zi ** 2 for ion in ions]
        Zeffx = np.sum([x * y for x, y in zip(ninormz, Ziz2)])
        
        for i, (key, value) in enumerate(self['meta'].items(), 1):
            print (value.__class__) 
            if not value.__class__ == list and not value.__class__ == np.ndarray:
                if i in [14, 15, 16, 17, 18, 19, 21, 22, 23, 24, 25, 26, 27, 28]:
                    print ('dup!')
                    bytes = array.array('d', [value]*len(scan_value))
                else:
                    bytes = array.array('d', [value])
            else:
                bytes = array.array('d', value)
            with open('input/p' + str(i) + '.bin', 'wb') as file:
                bytes.tofile(file)
#for i=1:scann  	   

#   Zeffx(i)=sum(ninormz(i,:).*Zi(i,:).^2); #calculate Zeff
#end
#
#figure; 
#disp('Figures displayed such that the input data can be verified by eye') 
#subplot(221)
#set(gca,'FontSize',18)
#plot(1:scann,Ati(:,1),'r',1:scann,Tix(:,1),'r--',1:scann,Ate,'b',1:scann,Tex,'b--','LineWidth',2)
#xlabel('Scan index')
#
#l1=legend('-R\nabla{T_i}/T_i','Ti','-R\nabla{T_e}/T_e','T_e');
#set(l1,'FontSize',12)
#legend boxoff
#grid on
#subplot(222)
#set(gca,'FontSize',18)
#plot(1:scann,Ane,'g',1:scann,Nex,'g--','LineWidth',2)
#xlabel('Scan index')
#
#l2=legend('-R\nabla{n_e}/n_e','n_e');
#set(l2,'FontSize',12)
#legend boxoff
#grid on
#subplot(223)
#set(gca,'FontSize',18)
#plot(1:scann,tau,'c',1:scann,Zeffx,'c--',1:scann,ft,'c-.','LineWidth',2)
#xlabel('Scan index')
#
#l3=legend('T_e/T_i','Z_{eff}','f_t');
#set(l3,'FontSize',12)
#legend boxoff
#grid on
#subplot(224)
#set(gca,'FontSize',18)
#plot(1:scann,smag,'m',1:scann,qx,'m--',1:scann,alphax,'m-.s','LineWidth',2)
#xlabel('Scan index')
#l4=legend('s','q','\alpha');
#set(l4,'FontSize',12)
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
sim = Sim()
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

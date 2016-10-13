"""
Copyright Dutch Institute for Fundamental Energy Research (2016)
Contributors: Karel van de Plassche (karelvandeplassche@gmail.com)
License: CeCILL v2.1
"""
import numpy as np
import matplotlib.pyplot as plt
from os import path

output_dir = 'output'
file_list = ['cftrans.dat', 'cke.dat', 'cki.dat', 'dfe_SI.dat', 'dfi_SI.dat',
             'ecoefs.dat', 'eef_cm.dat', 'eefETG_SI.dat', 'eef_GB.dat',
             'eef_SI.dat', 'epf_cm.dat', 'epfETG_SI.dat', 'epf_GB.dat',
             'epf_SI.dat', 'evf_cm.dat', 'evf_GB.dat', 'evf_SI.dat',
             'gam_GB.dat', 'gam_SI.dat', 'ief_cm.dat', 'ief_GB.dat',
             'ief_SI.dat', 'ipf_cm.dat', 'ipf_GB.dat', 'ipf_SI.dat',
             'ivf_cm.dat', 'ivf_GB.dat', 'ivf_SI.dat', 'modeflag.dat',
             'npol.dat', 'ome_GB.dat', 'ome_SI.dat', 'phi.dat', 'vce_SI.dat',
             'vci_SI.dat', 'vre_SI.dat', 'vri_SI.dat', 'vte_SI.dat',
             'vti_SI.dat']
file_list = [path.join(output_dir, file) for file in file_list]
             
def listify(file_path):
    with open(file_path) as file:
        array = []
        for line in file.readlines():
            words = line.split(None)
            line_array = []
            for word in words:
                line_array.append(float(word))
            array.append(line_array)
    return array
 
def plotScanWrapper(values, label=None, filter_zeros=True):
    if not np.all(values == 0.) or not filter_zeros:
        plt.plot(values, label=label)

def plotWrapper(data, plot_style):
    name, style, ylabel = plot_style
    print (data)
    print (name)
    print (style)
    plt.figure(name)
    plt.title(name)
    if style.startswith('w'):
        xlabel = r'wavenumber'
        for scan in data:
            plotScanWrapper(scan, label='temp')
    if style.startswith('i'):
        xlabel = r'scan value'
        data = data.T
        labels = ['D', 'Be', 'W']
        for label, scan in zip(labels, data):
            plotScanWrapper(scan, label, filter_zeros=False)
    if style.startswith('e'):
        xlabel = r'scan value'
        plotScanWrapper(data, filter_zeros=False)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()


plot_styles = {
    'gam_GB': ['Growth Rate [GB]', 'wavenumber', r'growth rate [$\sqrt{T_e/mi}/a$]'],
    'gam_SI': ['Growth Rate [SI]', 'wavenumber', r'growth rate [$\sqrt{T_e/mi}/a$]'],
    'ome_GB': ['Frequencies [GB]', 'wavenumber', r'frequency [?]'],
    'ome_SI': ['Frequencies [SI]', 'wavenumber', r'frequency [$s^-1$]'],
    
    'ief_GB': ['Ion Heat Conductivity [GB]', 'ions', r'heat conductivity [$\sqrt{mi}T_e^{1.5}/(q_e^2B^2a)$]'],
    'ief_SI': ['Ion Heat Flux [SI]', 'ions', r'heat flux [$W/m^2$]'],
    'eef_GB': ['Electron Heat Conductivity [GB]', 'electrons', r'heat conductivity [$\sqrt{mi}T_e^{1.5}/(q_e^2B^2a)$]'],
    'eef_SI': ['Electron Heat Flux [SI]', 'electrons', r'heat flux [$W/m^2$]'],
    'eefETG_SI': ['Electron Scale Heat Flux [SI]', 'electrons', r'heat flux [$W/m^2$]'],
    'ipf_GB': ['Ion Particle Diffusivity [GB]', 'ions', r'particle diffusifity [$\sqrt{mi}T_e^{1.5}/(q_e^2B^2a)$]'],
    'ipf_SI': ['Ion Particle Flux [SI]', 'ions', r'particle flux [$m^-2 s^-1$]'],
    'epf_GB': ['Electron Particle Diffusivity [GB]', 'electrons', r'particle diffusifity [$\sqrt{mi}T_e^{1.5}/(q_e^2B^2a)$]'],
    'epf_SI': ['Electron Particle Flux [GB]', 'electrons', r'particle flux [$m^-2 s^-1$]'],
    'epfETG_SI': ['Electron Scale Particle Flux [GB]', 'electrons', r'particle flux [$m^-2 s^-1$]'],
    'ivf_GB': ['Ion Momentum Diffusivity [GB]', 'ions', r'momentum diffusifity [$\sqrt{mi}T_e^{1.5}/(q_e^2B^2a)$]'],
    'ivf_SI': ['Ion Momentum Flux [SI]', 'ions', r'momentum flux [$N s m^-2 s^-1$]'],
    'evf_GB': ['Electron Momentum Diffusivity [GB]', 'electrons', r'momentum diffusifity [$\sqrt{mi}T_e^{1.5}/(q_e^2B^2a)$]'],
    'ivf_SI': ['Electron Momentum Flux [SI]', 'electrons', r'momentum flux [$N s m^-2 s^-1$]'],

}

prims = {
    'gam': ['growth rate', 'SI', 'GB'],
    'ome': ['frequencies', 'SI', 'GB'],
    'ef': ['heat', 'conductivity', 'flux'],
    'efETG': ['heat', 'conductivity', 'flux'],
    'pf': ['particle', 'diffusivity', 'flux'],
    'pfETG': ['particle', 'diffusivity', 'flux'],
    'vf': ['momentum', 'diffusivity', 'flux'],
    'df': ['', 'diffusivity', None],
    'vt': ['particle', 'thermopinch', None],
    'vr': ['particle', 'rotodiffusion pinch', None],
    'vc': ['particle', 'compressebility pinch', None],
    'chie': ['heat', 'conductivity', None],
    'ven': ['heat', 'thermopinch', None],
    'ver': ['heat', 'rotodiffusion pinch', None],
    'vec': ['heat', 'compressebility pinch', None]}
    

output = {}
for file_dir in file_list:
    __, temp = path.split(file_dir)
    name, __ = path.splitext(temp)
    base, __, unit = name.partition("_")
    if unit == "SI" and not (base.startswith("e") or base.startswith("i")) \
        and not (base == "gam" or base == "ome"):
        base = base[-1] + base[:-1]
    
    if base == "gam" or base == "ome":
        prim = base
        title = prims[prim][0] + " "
        if unit == "SI":
            title += prims[prim][1]
        elif unit == "GB":
            title += prims[prim][2]
        plotWrapper(np.array(listify(file_dir)), [title, 'w', "label"])
    elif not unit == "SI" and not unit == "GB":
        pass
    else:
        prim = base[1:]
        if base.startswith('i'):
            title = 'Ion '
        elif base.startswith('e'):
            title = 'Electron '
        title += prims[prim][0] + ' '
        if unit == "SI":
            title += prims[prim][1]
        elif unit == "GB":
            title += prims[prim][2]
        plotWrapper(np.array(listify(file_dir)), [title, base, "label"])
    output[name] = np.array(listify(file_dir))

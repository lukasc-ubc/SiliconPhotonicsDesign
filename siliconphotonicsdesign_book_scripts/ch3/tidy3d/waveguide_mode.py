# Waveguide mode simulations for strip waveguide, different widths
# using FlexCompute's tidy3d
# by Lukas Chrostowski, 2023

# tidy3d:
# https://docs.flexcompute.com/projects/tidy3d/en/stable/notebooks/ModeSolver.html

try:
    # install the packages in KLayout using SiEPIC-Tools function
    # Required packages for tidy3d:
    from SiEPIC.install import install
    install('plotly')
    install('pandas')
    install('packaging')
    install('defusedxml')
    install('numpy')
    install('toml')
    install('nest_asyncio')
    install('urllib3')

    # Import tidy3d, from SiEPIC-Tools package
    import SiEPIC.tidy3d.tidy3d as td
    from SiEPIC.tidy3d.tidy3d.constants import C_0
    from SiEPIC.tidy3d.tidy3d.plugins.mode.mode_solver import ModeSolver
    from SiEPIC.tidy3d.tidy3d import material_library

except:
    import tidy3d as td
    from tidy3d.constants import C_0
    from tidy3d.plugins.mode.mode_solver import ModeSolver
    from tidy3d import material_library

import plotly.express as px
import numpy as np


def wg_2D_draw(wg_width = 0.5, wg_height = 0.22):
    '''
    Setup a waveguide eigenmode simulation
    '''

    # waveguide materials
    medium_Si = material_library['cSi']['Li1993_293K']
    
    waveguide = td.Structure(
        geometry=td.Box(size=(wg_width, td.inf, wg_height)),
        medium=medium_Si,
    )
    
    return waveguide
    
def wg_2D(waveguide, wvl_um = 1.55, freq_span=0.2, num_freqs = 11, min_steps_per_wvl=40, Lx=2, Lz=1.5):
    '''
    Setup a 2D waveguide simulation
    '''

    # central frequency
    freq0 = C_0 / wvl_um
    fwidth = freq0 * freq_span
    
    # size of simulation domain
#    Lx, Ly, Lz = 2, 0, 1.5  # width, length, height
#    Lx, Ly, Lz = 4, 0, 3  # width, length, height
            
    # automatic grid specification
    grid_spec = td.GridSpec.auto(
        min_steps_per_wvl=min_steps_per_wvl, 
        wavelength=wvl_um,
    )

    sim = td.Simulation(
        size=(Lx, 0, Lz),
        grid_spec=grid_spec,
        structures=[waveguide],
        run_time = 1e-12,
        symmetry = (-1,0,0),  # for TE polarization
        medium = material_library['SiO2']['Horiba'],
        boundary_spec=td.BoundarySpec.all_sides(boundary=td.PECBoundary()),
    )
    
    plane = td.Box(center=(0, 0, 0), size=(4, 0, 3.5))
    
    mode_spec = td.ModeSpec(
        num_modes=1,
#        target_neff=2.0,
    )

    # frequency points
    if num_freqs>1:
        # f0_ind = num_freqs // 2
        freqs = np.linspace(freq0 - fwidth / 2, freq0 + fwidth / 2, num_freqs)
    else:
        freqs = [freq0]
    
    mode_solver = ModeSolver(
        simulation=sim,
        plane=plane,
        mode_spec=mode_spec,
        freqs=freqs,
    )
    return mode_solver

def wg_2D_neff(mode_solver, PLOT=True):
    '''
    Waveguide effective index at central wavelength
    '''

    # frequency points
    f0_ind = len(mode_solver.freqs)//2
    freqs = mode_solver.freqs

    # solve
    mode_data = mode_solver.solve()
    n_eff = mode_data.n_eff  # real part of the effective mode index
    
    if PLOT:
        # Plot mode profile, using plotly
        Ex= abs(mode_data.Ex.isel(mode_index=0, f=f0_ind))[:,0,:]
        fig = px.imshow(Ex.transpose(), title='Eigenmode profile, at %s nm' % (C_0/freqs[f0_ind]*1000) )
        fig.show()

    return n_eff

def wg_2D_convergenceTest(npoints = 5, PLOT=True):
    ''' Perform convergence tests to determine simulation span, mesh
    '''
    
    # width convergence
    aLx = np.linspace(0.7, 4, npoints)
    neffs = np.zeros(len(aLx))
    for i in range(len(aLx)):
        waveguide = wg_2D_draw()
        mode_solver = wg_2D(waveguide, wvl_um = 1.55, num_freqs = 1, Lx=aLx[i])
        neffs[i]=wg_2D_neff(mode_solver, PLOT=False)

    if PLOT:
        # Plot mode0 frequency dependance
        fig = px.line(x=aLx, y=neffs-neffs[-1], 
            title = 'convergence test',
            labels={'x':'Span: width', 'y':'Effective Index difference'}, markers=True)
        fig.show()

    # height convergence
    aLz = np.linspace(0.5, 3, npoints)
    neffs = np.zeros(len(aLz))
    for i in range(len(aLz)):
        waveguide = wg_2D_draw()
        mode_solver = wg_2D(waveguide, wvl_um = 1.55, num_freqs = 1, Lz=aLz[i])
        neffs[i]=wg_2D_neff(mode_solver, PLOT=False)

    if PLOT:
        # Plot mode0 frequency dependance
        fig = px.line(x=aLz, y=neffs-neffs[-1], 
            title = 'convergence test',
            labels={'x':'Span: height', 'y':'Effective Index difference'}, markers=True)
        fig.show()

    # mesh convergence
    min_steps_per_wvl = np.linspace(20, 100, npoints)
    neffs = np.zeros(len(min_steps_per_wvl))
    for i in range(len(min_steps_per_wvl)):
        waveguide = wg_2D_draw()
        mode_solver = wg_2D(waveguide, wvl_um = 1.55, num_freqs = 1, min_steps_per_wvl=min_steps_per_wvl[i])
        neffs[i]=wg_2D_neff(mode_solver, PLOT=False)

    if PLOT:
        # Plot mode0 frequency dependance
        fig = px.line(x=min_steps_per_wvl, y=neffs-neffs[-1], 
            title = 'convergence test',
            labels={'x':'Span: min_steps_per_wvl', 'y':'Effective Index difference'}, markers=True)
        fig.show()


def wg_2D_neff_sweep_wavelength(mode_solver, PLOT=True):
    '''
    Sweep the wavelength, find waveguide dispersion properties, plot
    '''

    # frequency points
    f0_ind = len(mode_solver.freqs)//2
    freqs = mode_solver.freqs
    if len(freqs) < 3:
        raise Exception('Need at least three frequency points')

    # Solve
    mode_data = mode_solver.solve()
    n_eff = mode_data.n_eff  # real part of the effective mode index


    # poly fit
    lambdas = C_0/n_eff.f.to_numpy()*1e-6 # wavelength in meters
    lambda0=float(lambdas[f0_ind])
    x_fit = lambdas - lambda0
    y_fit = n_eff[:,0].to_numpy()
    if 0:
        print(lambda0)
        print(lambdas)
        print(x_fit)
        print(y_fit)
    a_fit = np.polyfit(x_fit, y_fit, 2)
    #print (a_fit)
    print("Fit function for effective index: neff(lambda) = %s + %s*(lambda-%s) + %s*(lambda-%s)^2" % 
        (a_fit[-1], a_fit[-2], lambda0, a_fit[-3], lambda0) );

    # convert to group index and dispersion
    neff = a_fit[-1]
    print("Effective index at %s: %s" % (lambda0,neff));
    ng = a_fit[-1] - a_fit[-2] * lambda0;
    print("Group index at %s: %s " %(lambda0, ng));
    D = -lambda0 / C_0 * 2 * a_fit[-3];
    print( "Dispersion [s/m^2] at %s: %s" %(lambda0,D));
    wg_dispersion = [neff, ng, D]


    if PLOT:
        # Plot mode0 frequency dependance
        fig = px.line(x=n_eff.f, y=n_eff[:,0], labels={'x':'Frequency', 'y':'Effective Index'}, markers=True)
        fig.show()

    return wg_dispersion

def wg_2D_neff_sweep_width(min_steps_per_wvl=20, PLOT=1):
    '''
    Sweep the width of the waveguide, and 
    extract waveguide dispersion versus width
    '''
    
#    widths = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
    widths = np.linspace(0.2, 0.7, 20)
    widths = np.linspace(0.4, 0.5, 20)
    
    neffs = np.zeros(len(widths))
    ngs = np.zeros(len(widths))
    disps = np.zeros(len(widths))
    
    for i in range(len(widths)):
        waveguide = wg_2D_draw(wg_width=widths[i])
        mode_solver = wg_2D(waveguide, 
            wvl_um = 1.55, 
            num_freqs = 3,
            freq_span=0.1,
            min_steps_per_wvl=min_steps_per_wvl,
        )
        neffs[i],ngs[i],disps[i] = wg_2D_neff_sweep_wavelength(mode_solver, PLOT=False)

    if PLOT:
        # Plot mode0 frequency dependance
        fig = px.line(x=widths, y=neffs, labels={'x':'Widths', 'y':'Effective Index'}, markers=True)
        fig.show()
        fig = px.line(x=widths, y=ngs, labels={'x':'Widths', 'y':'Group Index'}, markers=True)
        fig.show()
        fig = px.line(x=widths, y=disps, labels={'x':'Widths', 'y':'Dispersion'}, markers=True)
        fig.show()


# Example simulations

run = [0,3]

if 0 in run:
    # convergence test on neff
    wg_2D_convergenceTest(npoints = 30)
    
if 1 in run:
    # plot mode profile
    waveguide = wg_2D_draw()
    mode_solver = wg_2D(waveguide, wvl_um = 1.55, min_steps_per_wvl=40, num_freqs = 1)
    wg_2D_neff(mode_solver)

if 2 in run:
    # get waveguide dispersion
    waveguide = wg_2D_draw()
    mode_solver = wg_2D(waveguide, wvl_um = 1.55, num_freqs = 3)
    wg_dispersion = wg_2D_neff_sweep_wavelength(mode_solver, PLOT=True)

if 3 in run:
    # waveguide width sweep
    wg_2D_neff_sweep_width(min_steps_per_wvl=100)



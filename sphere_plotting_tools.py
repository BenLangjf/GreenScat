# Sphere Plotting tools
# Ben Lang 2025

import numpy as np
import matplotlib.pyplot as plt
import cmocean # ---> My old favourite for colourmaps
import matplotlib.cm as cm
import matplotlib as mpl


def theta_phi_to_E( theta, phi ):
    # Turns theta/phi coordinates on the Poincare sphere into an E-field.
    
    E = [ np.cos(theta/2) * np.cos(phi/2) - 1j * np.sin(theta/2) * np.sin(phi/2),
         np.sin(theta/2) * np.cos(phi/2) + 1j * np.cos(theta/2) * np.sin(phi/2)]
    
    return E
    


def make_sphere_plot( n_theta, n_phi, function_to_apply ):
    
    theta = np.linspace(0, 2*np.pi, n_theta )
    phi   = np.linspace( -np.pi/2,   np.pi/2, n_phi )
    
    
    p_grid, t_grid = np.meshgrid( phi, theta)
    
    data = function_to_apply( t_grid, p_grid )
    
    # Get Streams of data phase angle gradient
    dp, dt = np.gradient( np.angle(data) )
    
    fi, ax = plt.subplots(figsize = (5,5))
    stream = ax.streamplot( p_grid, t_grid, dt, dp )
    paths = stream.lines.get_paths()
    segments = stream.lines.get_segments()

    
    # Get XYZ coordinates
    def get_xyz(theta, phi):
        return [ np.cos(theta)*np.cos(phi), np.sin(theta)*np.cos(phi), np.sin(phi) ]
    
    x, y, z = get_xyz( t_grid, p_grid )
    
    # Transform data into colours
    norm = mpl.colors.Normalize(vmin=-np.pi, vmax=+np.pi)
    mapper = cm.ScalarMappable( norm = norm, cmap= cmocean.cm.phase )
    cols = mapper.to_rgba( np.angle(data) )
    
    # Use Amplitude to get Alpha Channel
    # Bring between 0 and 1:
    amplitudes = np.abs(data)
    cols[:,:,3] = amplitudes / amplitudes.max()
    
    fig, ax = plt.subplots(figsize = (5,5))
    ax.imshow( cols )
    
    # However, matplotlib cannot do alpha channel in 3D.
    # So I need to use that channel to "whiten" everything.
    # When amplitude is 0, cols need to go to white (1,1,1)
    
    white = np.ones( cols[:,:,0].shape )
    
    cols[:,:,0] = cols[:,:,0] * cols[:,:,3] + white * (1 - cols[:,:,3])
    cols[:,:,1] = cols[:,:,1] * cols[:,:,3] + white * (1 - cols[:,:,3])
    cols[:,:,2] = cols[:,:,2] * cols[:,:,3] + white * (1 - cols[:,:,3])
    
    fig, ax = plt.subplots(figsize = (5,5))
    ax.imshow( cols[:,:,:3] )
    
    fig, ax3d = plt.subplots(subplot_kw={"projection": "3d"})
    ax3d.mouserotationstyle: trackball
    plt.rcParams['patch.edgecolor'] = 'none'
    ax3d.plot_surface(x, y, z, facecolors = cols[:,:,:3], shade = False, linewidth = 0, antialiased = True )
    ax3d.set_xlim3d([-1.1, 1.1])
    ax3d.set_ylim3d([-1.1, 1.1])
    ax3d.set_zlim3d([-1.1, 1.1])
    ax3d.set_axis_off
    
    
    # # Get XYZ coordinates of paths
    paths_xyz = [ get_xyz(path.vertices[:,1], path.vertices[:,0]) for path in paths ]
    #paths_xyz = np.array(paths_xyz)
    
    fig, ax3d = plt.subplots(subplot_kw={"projection": "3d"})
    ax3d.mouserotationstyle: trackball
    plt.rcParams['patch.edgecolor'] = 'none'
    for path_xyz in paths_xyz:
        ax3d.plot( *path_xyz, 'k' )
        
    # Urghh, in MATLAB this "just works". In python you have to jump through 100 hoops
    #ax3d.plot_surface(0.9*x, 0.9*y, 0.9*z, facecolor = 'w', shade = False, linewidth = 0, antialiased = True, alpha = 0.9 )
    ax3d.set_xlim3d([-1.1, 1.1])
    ax3d.set_ylim3d([-1.1, 1.1])
    ax3d.set_zlim3d([-1.1, 1.1])
    ax3d.set_axis_off
    

from MGS import *


@np.vectorize
def scattering_from_various_dipoles_C_point( theta, phi ):
    # Use Vecotrise so that when fed the meshgrids it does them element-wise. (map)

    vg = 0.1
    mag = (0.5 / vg)

    # Isotropically Polarisable dipole
    S = multilevel_system([10], [0])
    S.add_transition( 0, 0,  theta_phi_to_E(theta, phi) )
    
    G = WG_Greens_function( [[1, 1j]], [mag], ["WG"])
    G.add_mode( [1, 0], 0.3, False, "H loss" )
    G.add_mode( [0, 1], 0.3, False, "V loss" )
    
    scat, space = long_photon_scattering(G, S, [1], [1, 0, 0, 0], freq=10)
    
    return scat[1]


@np.vectorize
def scattering_low_loss_isotropic( theta, phi ):
    # Use Vecotrise so that when fed the meshgrids it does them element-wise. (map)

    vg = 0.1
    mag = (0.5 / vg)

    # Isotropically Polarisable dipole
    S = multilevel_system([10, 10], [0])
    #S.add_transition( 0, 0, [1, +1j] )
    #S.add_transition( 1, 0, [1, -1j] )
    S.add_transition( 0, 0, [1, 0] )
    S.add_transition( 1, 0, [0, 1] )
    
    G = WG_Greens_function( [theta_phi_to_E(theta, phi)], [mag], ["WG"])
    G.add_mode( [1, 0], 0.01, False, "H loss" )
    G.add_mode( [0, 1], 0.01, False, "V loss" )
    
    scat, space = long_photon_scattering(G, S, [1], [1, 0, 0, 0], freq=10)
    
    return scat[1]




#make_sphere_plot( 200, 100, scattering_from_various_dipoles_C_point )
make_sphere_plot( 200, 100, scattering_low_loss_isotropic )
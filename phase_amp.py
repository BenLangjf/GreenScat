import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage
import matplotlib as mpl

def phase_amp( field, xvec, yvec = 'same', ax = 'none', shift = (0,0), half_off = False, root = False, aspect = 1 ):
    # Phase and Amplitude plot of complex field.
    
    if yvec =='same':
        yvec = xvec
    
    if ax == 'none':
        fig, ax = plt.subplots(figsize=(7,7))
    
    x, y = shift
    
    phases = np.angle( field )
    amplitudes = np.abs( field )
    
    import cmocean # ---> My old favourite for colourmaps
    import matplotlib.cm as cm
    
    
    if half_off:
        # Alter the Colour-map on one half of the area to highlight small stuff.
        sqrt_amps = amplitudes ** (1/3)
        sqrt_amps *= ( amplitudes.max() / sqrt_amps.max() )
        
        for n in range(len(amplitudes)):
            amplitudes[:n,n] = sqrt_amps[:n,n]
    elif root:
        amplitudes = np.sqrt( amplitudes )
        
    # Use phase to set colour
    norm = mpl.colors.Normalize(vmin=-np.pi, vmax=+np.pi)
    mapper = cm.ScalarMappable( norm = norm, cmap= cmocean.cm.phase )
    #mapper = cm.ScalarMappable( norm = lambda x: (x+np.pi) / (2*np.pi), cmap= cmocean.cm.phase )
    cols = mapper.to_rgba( phases )
    
    # Use Amplitude to get Alpha Channel
    # Bring between 0 and 1:
    amplitudes = amplitudes / amplitudes.max()
    cols[:,:,3] = amplitudes
    
    ax.imshow( cols, extent = [xvec[0] + x, xvec[-1] + x, yvec[0] + y, yvec[-1] + y ], origin = 'lower', zorder = 0, aspect = aspect)
    #ax.axis('equal')
    ax.set_xlim( [xvec[0], xvec[-1]] )
    ax.set_ylim( [yvec[0], yvec[-1]] )
    
    if half_off:
        ax.plot( [ xvec[0], xvec[-1]], [yvec[0], yvec[-1]], 'k-' , zorder = 1)
        
    
    return ax
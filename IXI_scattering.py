
from MGS import *

####### Scattering from an Isotropically polarisable dipole #######

vg = 0.1
mag = (0.5 / vg)

# Try IXI system (Vogit QD)
S = multilevel_system([10, 10], [0, 0])
S.add_transition( 0, 0, [1, 0] )
S.add_transition( 1, 1, [1, 0] )
S.add_transition( 0, 1, [0, 1j] )
S.add_transition( 1, 0, [0, 1j] )

thetas = np.linspace(0, np.pi, 300)
data = []

mats = []

for theta in thetas:
    
    G = WG_Greens_function( [[np.cos(theta), 1j*np.sin(theta)]], [mag], ["WG"])
    G.add_mode( [1, 0], 0.2, False, "H loss" )
    G.add_mode( [0, 1], 0.2, False, "V loss" )
    
    
    scat, space = long_photon_scattering(G, S, [1, 0], [1, 0, 0, 0], freq=10)
    data.append( scat )
    
    
data = np.array(data)
data = data[:,:,0]


fig, ax = plt.subplots(figsize = (5,5))
ax.plot( thetas, data[:,0], 'b-', label = space[0] )
ax.plot( thetas, (data[:,1]), 'b--' , label = space[1])
ax.plot( thetas, (data[:,2]), 'r-' , label = space[2])
ax.plot( thetas, (data[:,3]), 'r--' , label = space[3])

# Note : I am taking the abs of the reflection, as their is an 
# arbitary transformation on the phase. (I can multiply the Green's function
# polarisation by exp(i theta), which changes the reflection phase by a factor
# exp(2i theta), but the whole transformation has no measureable effect, and is 
# essentially just a basis change.)
# Transmission is different - its phase is meaningful, as it is measured
# relative to the phase with no scatterer, instead of an arbitary benchmark.

# Although I don't have to do this. Saying the polarisation is [0, 1j] does imply
# that the scattering atom is a quater wave further away from your arbitary reference
# point that saying the polarisation is [0, 1].

for theta in thetas[::15]:
    draw_ellipse( np.cos(theta), 1j*np.sin(theta), centre = [theta, -1.1], ax = ax, scale=0.12, norm=True, color = 'r' )


ax.set_xticks([0, np.pi], ['0', r'$\pi$'])
ax.set_xlabel( r'$\theta$' )
ax.set_ylabel('Amplitude')

ax.set_xlim([0, np.pi])

ax.legend()

fig.savefig("IXI_Scattering.pdf", bbox_inches='tight')

S.represent()


from MGS import *

####### Scattering from an Isotropically polarisable dipole #######

vg = 0.1
mag = (0.5 / vg)

# Try Lambda system
S = multilevel_system([10], [0, 0])
S.add_transition( 0, 0, [1, +1j] )
S.add_transition( 0, 1, [1, -1j] )

S2 = multilevel_system([10], [0, 0])
S2.add_transition( 0, 0, [1, +2j] )
S2.add_transition( 0, 1, [1, -2j] )

thetas = np.linspace(0, np.pi, 300)

data = []
data2 = []

for theta in thetas:
    
    G = WG_Greens_function( [[np.cos(theta), 1j*np.sin(theta)]], [mag], ["WG"])
    G.add_mode( [1, 0], 0.2, False, "H loss" )
    G.add_mode( [0, 1], 0.2, False, "V loss" )
    
    
    scat, space = long_photon_scattering(G, S, [1, 0], [1, 0, 0, 0], freq=10)
    data.append( scat )
    
    scat, space = long_photon_scattering(G, S2, [1, 0], [1, 0, 0, 0], freq=10)
    data2.append( scat )
    
    
data = np.array(data)
data = data[:,:,0]

data2 = np.array(data2)
data2 = data2[:,:,0]

fig, axes = plt.subplots(2, 1, figsize = (5,10))
axes[0].plot( thetas, data[:,0], 'b-', label = space[0] )
axes[0].plot( thetas, (data[:,1]), 'b--' , label = space[1])
axes[0].plot( thetas, (data[:,2]), 'r-' , label = space[2])
axes[0].plot( thetas, (data[:,3]), 'r--' , label = space[3])


for theta in thetas[::15]:
    draw_ellipse( np.cos(theta), 1j*np.sin(theta), centre = [theta, -1.1], ax = axes[0], scale=0.12, norm=True, color = 'r' )


axes[0].set_xticks([0, np.pi], ['0', r'$\pi$'])
axes[0].set_xlabel( r'$\theta$' )
axes[0].set_ylabel('Amplitude')

axes[0].set_xlim([0, np.pi])

axes[0].legend()


axes[1].plot( thetas, data2[:,0], 'b-', label = space[0] )
axes[1].plot( thetas, (data2[:,1]), 'b--' , label = space[1])
axes[1].plot( thetas, (data2[:,2]), 'r-' , label = space[2])
axes[1].plot( thetas, (data2[:,3]), 'r--' , label = space[3])

for theta in thetas[::15]:
    draw_ellipse( np.cos(theta), 1j*np.sin(theta), centre = [theta, -1.1], ax = axes[1], scale=0.12, norm=True, color = 'r' )


axes[1].set_xticks([0, np.pi], ['0', r'$\pi$'])
axes[1].set_xlabel( r'$\theta$' )
axes[1].set_ylabel('Amplitude')

axes[1].set_xlim([0, np.pi])

axes[1].legend()


fig.savefig("Lambda_Scattering.pdf", bbox_inches='tight')

S.represent()
S2.represent()
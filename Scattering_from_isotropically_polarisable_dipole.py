
from MGS import *

####### Scattering from an Isotropically polarisable dipole #######

vg = 0.1
mag = (0.5 / vg)

# Isotropically Polarisable dipole
S = multilevel_system([10, 10], [0])
S.add_transition( 0, 0, [1, +1j] )
S.add_transition( 1, 0, [1, -1j] )

thetas = np.linspace(0, np.pi, 300)
data1 = []
data2 = []

mats = []

for theta in thetas:
    
    G = WG_Greens_function( [[np.cos(theta), 1j*np.sin(theta)]], [mag], ["WG"])
    G.add_mode( [1, 0], 0.2, False, "H loss" )
    G.add_mode( [0, 1], 0.2, False, "V loss" )
    
    G2 = WG_Greens_function( [[np.cos(theta), 1j*np.sin(theta)]], [mag], ["WG"])
    G2.add_mode( [1, 0], 0.003, False, "H loss" )
    G2.add_mode( [0, 1], 0.003, False, "V loss" )
    
    scat, space = long_photon_scattering(G, S, [1], [1, 0, 0, 0], freq=10)
    data1.append( scat )
    
    scat, space = long_photon_scattering(G2, S, [1], [1, 0, 0, 0], freq=10)
    data2.append( scat )
    
    # Just Gamma matrix, before inversion stuff.
    mats.append( S.inner( G2 ) )
    

data1 = np.array(data1)
data1 = data1[:,:,0]

data2 = np.array(data2)
data2 = data2[:,:,0]

fig, ax = plt.subplots(figsize = (5,5))
ax.plot( thetas, data1[:,0], 'b-', label = r"$t$" )
ax.plot( thetas, np.abs(data1[:,1]), 'r-' , label = r'$r$')
ax.plot( thetas, data2[:,0], 'b--' , label = r"$t$ low loss" )
ax.plot( thetas, np.abs(data2[:,1]), 'r--' , label = r"$r$ low loss" )

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

# Draw polarisation ellipse of "effective dipole" (read out from loss modes)
for theta, data2_point in zip(thetas[::15], data2[::15]):
    draw_ellipse( data2_point[2], data2_point[3], centre = [theta, -1.3], ax = ax, scale=0.12, norm=True, color = 'k' )



ax.set_xticks([0, np.pi], ['0', r'$\pi$'])
ax.set_xlabel( r'$\theta$' )
ax.set_ylabel('Amplitude')

ax.set_xlim([0, np.pi])

ax.legend()

fig.savefig("Scattering_Isotropic.pdf", bbox_inches='tight')


#%% Check something. Is it just like a two level system but one that has dipole matching the polarisation.

data_I = []

for theta in thetas:
    
    S_I = multilevel_system([10], [0])
    S_I.add_transition( 0, 0, [np.sin(theta), 1j*np.cos(theta)] )
    
    G = WG_Greens_function( [[np.cos(theta), 1j*np.sin(theta)]], [mag], ["WG"])
    G.add_mode( [1, 0], 0.003, False, "H loss" )
    G.add_mode( [0, 1], 0.003, False, "V loss" )
    
    scat, space = long_photon_scattering(G, S_I, [1], [1, 0, 0, 0], freq=10)
    data_I.append( scat )

    
data_I = np.array(data_I)
data_I = data_I[:,:,0]

ax.plot( thetas, data_I[:,0], 'bx', label = r"$t$" )
ax.plot( thetas, data_I[:,1], 'rx' , label = r'$r$')


#%% That is only one line going from North pole to South on our Poincare sphere.
# We need to look at some different lines as well!


thetas = np.linspace(0, np.pi, 300)
data1 = []
data2 = []

mats = []

for theta in thetas:
    
    G = WG_Greens_function( [[np.cos(theta), (1/np.sqrt(2))*(1+1j)*np.sin(theta)]], [mag], ["WG"])
    G.add_mode( [1, 0], 0.2, False, "H loss" )
    G.add_mode( [0, 1], 0.2, False, "V loss" )
    
    G2 = WG_Greens_function( [[np.cos(theta), (1/np.sqrt(2))*(1+1j)*np.sin(theta)]], [mag], ["WG"])
    G2.add_mode( [1, 0], 0.003, False, "H loss" )
    G2.add_mode( [0, 1], 0.003, False, "V loss" )
    
    scat, space = long_photon_scattering(G, S, [1], [1, 0, 0, 0], freq=10)
    data1.append( scat )
    
    scat, space = long_photon_scattering(G2, S, [1], [1, 0, 0, 0], freq=10)
    data2.append( scat )
    
    # Just Gamma matrix, before inversion stuff.
    mats.append( S.inner( G2 ) )
    

data1 = np.array(data1)
data1 = data1[:,:,0]

data2 = np.array(data2)
data2 = data2[:,:,0]

fig, ax = plt.subplots(figsize = (5,5))
ax.plot( thetas, data1[:,0], 'b-', label = r"$t$" )
ax.plot( thetas, np.abs(data1[:,1]), 'r-' , label = r'$r$')
ax.plot( thetas, data2[:,0], 'b--' , label = r"$t$ low loss" )
ax.plot( thetas, np.abs(data2[:,1]), 'r--' , label = r"$r$ low loss" )

# Note : I am taking the abs of the reflection, as their is an 
# arbitary transformation on the phase. (I can multiply the Green's function
# polarisation by exp(i theta), which changes the reflection phase by a factor
# exp(2i theta), but the whole transformation has no measureable effect, and is 
# essentially just a basis change.)
# Transmission is different - its phase is meaningful, as it is measured
# relative to the phase with no scatterer, instead of an arbitary benchmark.

for theta in thetas[::15]:
    draw_ellipse( np.cos(theta), (1/np.sqrt(2))*(1+1j)*np.sin(theta), centre = [theta, -1.1], ax = ax, scale=0.12, norm=True, color = 'r' )

# Draw polarisation ellipse of "effective dipole" (read out from loss modes)
for theta, data2_point in zip(thetas[::15], data2[::15]):
    draw_ellipse( data2_point[2], data2_point[3], centre = [theta, -1.3], ax = ax, scale=0.12, norm=True, color = 'k' )



ax.set_xticks([0, np.pi], ['0', r'$\pi$'])
ax.set_xlabel( r'$\theta$' )
ax.set_ylabel('Amplitude')

ax.set_xlim([0, np.pi])
ax.text( 1, 1, 'But only real parts plotted! Important data unseen.' )

ax.legend()



#%% That is only one line going from North pole to South on our Poincare sphere.
# We need to look at some different lines as well!
# The Equator is presumably boring?


thetas = np.linspace(0, np.pi, 300)
data1 = []
data2 = []

mats = []

for theta in thetas:
    
    G = WG_Greens_function( [[np.cos(theta), np.sin(theta)]], [mag], ["WG"])
    G.add_mode( [1, 0], 0.2, False, "H loss" )
    G.add_mode( [0, 1], 0.2, False, "V loss" )
    
    G2 = WG_Greens_function( [[np.cos(theta), np.sin(theta)]], [mag], ["WG"])
    G2.add_mode( [1, 0], 0.003, False, "H loss" )
    G2.add_mode( [0, 1], 0.003, False, "V loss" )
    
    scat, space = long_photon_scattering(G, S, [1], [1, 0, 0, 0], freq=10)
    data1.append( scat )
    
    scat, space = long_photon_scattering(G2, S, [1], [1, 0, 0, 0], freq=10)
    data2.append( scat )
    
    # Just Gamma matrix, before inversion stuff.
    mats.append( S.inner( G2 ) )
    

data1 = np.array(data1)
data1 = data1[:,:,0]

data2 = np.array(data2)
data2 = data2[:,:,0]

fig, ax = plt.subplots(figsize = (5,5))
ax.plot( thetas, data1[:,0], 'b-', label = r"$t$" )
ax.plot( thetas, data1[:,1], 'r-' , label = r'$r$')
ax.plot( thetas, data2[:,0], 'b--' , label = r"$t$ low loss" )
ax.plot( thetas, data2[:,1], 'r--' , label = r"$r$ low loss" )

for theta in thetas[::15]:
    draw_ellipse( np.cos(theta), np.sin(theta), centre = [theta, -1.1], ax = ax, scale=0.12, norm=True, color = 'r' )

# Draw polarisation ellipse of "effective dipole" (read out from loss modes)
for theta, data2_point in zip(thetas[::15], data2[::15]):
    draw_ellipse( data2_point[2], data2_point[3], centre = [theta, -1.3], ax = ax, scale=0.12, norm=True, color = 'k' )



ax.set_xticks([0, np.pi], ['0', r'$\pi$'])
ax.set_xlabel( r'$\theta$' )
ax.set_ylabel('Amplitude')

ax.set_xlim([0, np.pi])

ax.legend()



#%% That is only one line going from North pole to South on our Poincare sphere.
# We need to look at some different lines as well!


thetas = np.linspace(0, np.pi, 300)
data1 = []
data2 = []

mats = []

for phi in thetas:
    
    G = WG_Greens_function( [[(1/np.sqrt(2)), (1/np.sqrt(2)) * np.exp(1j*phi)]], [mag], ["WG"])
    G.add_mode( [1, 0], 0.2, False, "H loss" )
    G.add_mode( [0, 1], 0.2, False, "V loss" )
    
    G2 = WG_Greens_function( [[(1/np.sqrt(2)), (1/np.sqrt(2)) * np.exp(1j*phi)]], [mag], ["WG"])
    G2.add_mode( [1, 0], 0.003, False, "H loss" )
    G2.add_mode( [0, 1], 0.003, False, "V loss" )
    
    scat, space = long_photon_scattering(G, S, [1], [1, 0, 0, 0], freq=10)
    data1.append( scat )
    
    scat, space = long_photon_scattering(G2, S, [1], [1, 0, 0, 0], freq=10)
    data2.append( scat )
    
    # Just Gamma matrix, before inversion stuff.
    mats.append( S.inner( G2 ) )
    

data1 = np.array(data1)
data1 = data1[:,:,0]

data2 = np.array(data2)
data2 = data2[:,:,0]

fig, ax = plt.subplots(figsize = (5,5))
ax.plot( thetas, data1[:,0], 'b-', label = r"$t$" )
ax.plot( thetas, data1[:,1], 'r-' , label = r'$r$')
ax.plot( thetas, data2[:,0], 'b--' , label = r"$t$ low loss" )
ax.plot( thetas, data2[:,1], 'r--' , label = r"$r$ low loss" )

for phi in thetas[::15]:
    draw_ellipse( 1/np.sqrt(2), 1/np.sqrt(2)*np.exp(1j*phi), centre = [phi, -1.1], ax = ax, scale=0.12, norm=True, color = 'r' )

# Draw polarisation ellipse of "effective dipole" (read out from loss modes)
for theta, data2_point in zip(thetas[::15], data2[::15]):
    draw_ellipse( data2_point[2], data2_point[3], centre = [theta, -1.3], ax = ax, scale=0.12, norm=True, color = 'k' )



ax.set_xticks([0, np.pi], ['0', r'$\pi$'])
ax.set_xlabel( r'$\theta$' )
ax.set_ylabel('Amplitude')

ax.set_xlim([0, np.pi])

ax.legend()


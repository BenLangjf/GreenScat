
from MGS import *

####### TESTS #######


vg = 0.1
mag = (0.5 / vg)

#%% Linear Dipole Two-level system Reflection

G = WG_Greens_function( [[1, 1j]], [mag], ["WG"])
G.add_mode( [1, 0], 0.1, False, "H loss" )
G.add_mode( [0, 1], 0.1, False, "V loss" )


S = multilevel_system([10], [0])
S.add_transition( 0, 0, [1, 0] )

scat, space = long_photon_scattering(G, S, [1], [1, 0, 0, 0], freq=10)

print('Linear dipole reflection')
for r, name in zip([complex(i[0,0]) for i in scat], space):
    print( name, " " * (24-len(name)), r)
    
print("Norm = ", np.linalg.norm( np.array(scat) ), "\n" )
## T1 ended

#%% Circular Dipole Two-level system with circular polarisation, transmission.

S = multilevel_system([10], [0])
S.add_transition( 0, 0, [1, +1j] )

scat, space = long_photon_scattering(G, S, [1], [1, 0, 0, 0], freq=10)

print('Circular dipole phase shift')
for r, name in zip([complex(i[0,0]) for i in scat], space):
    print( name, " " * (24-len(name)), r)
    
print("Norm = ", np.linalg.norm( np.array(scat) ), "\n" )
## T2 ended

scat, space = long_photon_scattering(G, S, [1], [0, 1, 0, 0], freq=10)

print('Circular dipole backwards')
for r, name in zip([complex(i[0,0]) for i in scat], space):
    print( name, " " * (24-len(name)), r)
    
print("Norm = ", np.linalg.norm( np.array(scat) ), "\n" )
## T3 ended

#%% lambda system, circle dipole c point.
# Results in "toggle swtitch" (see doi:10.1126/science.1254699)

S = multilevel_system([10], [0, 0])
S.add_transition( 0, 0, [1, +1j] )
S.add_transition( 0, 1, [1, -1j] )

scat, space = long_photon_scattering(G, S, [1, 0], [1, 0, 0, 0], freq=10)

print('Lambda System, circular dipoles, c point. Should toggle-reflect.')
for r, name in zip([complex(i[0,0]) for i in scat], space):
    print( name, " " * (24-len(name)), r)
    
print("Norm = ", np.linalg.norm( np.array(scat) ), "\n" )
## T4 ended

#%%  II level system (charged QD)
S = multilevel_system([10, 10], [0, 0])
S.add_transition( 0, 0, [1, +1j] )
S.add_transition( 1, 1, [1, -1j] )

# Results should match 2 level system above:
scat, space = long_photon_scattering(G, S, [1, 0], [1, 0, 0, 0], freq=10)

print('II System, circular dipoles, should match circular dipole phase shift.')
for r, name in zip([complex(i[0,0]) for i in scat], space):
    print( name, " " * (24-len(name)), r)
    
print("Norm = ", np.linalg.norm( np.array(scat) ), "\n" )
## T5 ended

#%%  V level system
S = multilevel_system([10, 10], [0])
S.add_transition( 0, 0, [1, +1j] )
S.add_transition( 1, 0, [1, -1j] )


scat, space = long_photon_scattering(G, S, [1], [1, 0, 0, 0], freq=10)

print('V system phase shift')
for r, name in zip([complex(i[0,0]) for i in scat], space):
    print( name, " " * (24-len(name)), r)
    
print("Norm = ", np.linalg.norm( np.array(scat) ), "\n" )
## T6 ended


#%% linear V level system
S = multilevel_system([10, 10], [0])
S.add_transition( 0, 0, [1, 0] )
S.add_transition( 1, 0, [0, 1] )


scat, space = long_photon_scattering(G, S, [1], [1, 0, 0, 0], freq=10)

print('Linear V system phase shift')
for r, name in zip([complex(i[0,0]) for i in scat], space):
    print( name, " " * (24-len(name)), r)
    
print("Norm = ", np.linalg.norm( np.array(scat) ), "\n" )
## T7 ended

#%% V-level at elliptical point
vg = 0.1
mag = (0.5 / vg)

G = WG_Greens_function( [[2, 1j]], [mag], ["WG"])
G.add_mode( [1, 0], 0.1, False, "H loss" )
G.add_mode( [0, 1], 0.1, False, "V loss" )

scat, space = long_photon_scattering(G, S, [1], [1, 0, 0, 0], freq=10)

print('V level elliptical point')
for r, name in zip([complex(i[0,0]) for i in scat], space):
    print( name, " " * (24-len(name)), r)

# Note that the "effective dipole" is almost the cross-ellipse configuration
# For transmission with phase shift.

print("Norm = ", np.linalg.norm( np.array(scat) ), "\n" )
## T8 ended

#%% Test detuning, reflection
vg = 0.1
mag = (0.5 / vg)

G = WG_Greens_function( [[1, 0]], [mag], ["WG"])
G.add_mode( [1, 0], 0.1, False, "H loss" )
G.add_mode( [0, 1], 0.1, False, "V loss" )

S = multilevel_system([1000], [0])
S.add_transition( 0, 0, [1, 0] )

freqs = np.linspace(900,1100, 101)
dat = []
for freq in freqs:
    scat, space = long_photon_scattering(G, S, [1], [1, 0, 0, 0], freq=freq)
    dat.append( complex(scat[1,0]) )

dat = np.array(dat)

fig, ax = plt.subplots(figsize = (5,5))
ax.plot( freqs, np.abs(dat) )
ax2 = ax.twinx()
ax2.plot( freqs, np.angle(dat), 'r', label = "phase")
## T9 ended

#%% Test detuning, transmission
G = WG_Greens_function( [[1, 1j]], [mag], ["WG"])
G.add_mode( [1, 0], 0.1, False, "H loss" )
G.add_mode( [0, 1], 0.1, False, "V loss" )

S = multilevel_system([1000], [0])
S.add_transition( 0, 0, [1, 1j] )

freqs = np.linspace(900,1100, 101)
dat = []
for freq in freqs:
    scat, space = long_photon_scattering(G, S, [1], [1, 0, 0, 0], freq=freq)
    dat.append( complex(scat[0,0]) )

dat = np.array(dat)

fig, ax = plt.subplots(figsize = (5,5))
ax.plot( freqs, np.abs(dat) )

ax2 = ax.twinx()
ax2.plot( freqs, np.angle(dat), 'r' , label = "phase")


#%% Various Ellipses
# To match fig.3 of "Perfect Chirality with imperfect polarisation" 
# https://doi.org/10.1103/PhysRevLett.128.073602

vg = 0.1
mag = (0.5 / vg)

# Circular Polarisation
G = WG_Greens_function( [[1, 1j]], [mag], ["WG"])
G.add_mode( [1, 0], 0.05, False, "H loss" )
G.add_mode( [0, 1], 0.05, False, "V loss" )

# Elliptical Polarisation
# Note, the paper contains a typo, where it says the (normalised) polarsation is
# (13, 1j) / root(10) it shouls day (3, 1j)/root(10).

G2 = WG_Greens_function( [[1, 3j]], [mag], ["WG"])
G2.add_mode( [1, 0], 0.05, False, "H loss" )
G2.add_mode( [0, 1], 0.05, False, "V loss" )
    

thetas = np.linspace(0, np.pi, 300)
data1 = []
data2 = []

for theta in thetas:
    
    # Polarisable dipole
    S = multilevel_system([10, 10], [0])
    S.add_transition( 0, 0, [ np.cos(theta), 1j*np.sin(theta) ] )
    
    scat, space = long_photon_scattering(G, S, [1], [1, 0, 0, 0], freq=10)
    data1.append( scat )
    
    scat, space = long_photon_scattering(G2, S, [1], [1, 0, 0, 0], freq=10)
    data2.append( scat )
    

data1 = np.array(data1)
data1 = data1[:,:,0]

data2 = np.array(data2)
data2 = data2[:,:,0]

fig, axes = plt.subplots(2, 1, figsize = (5,10))
axes[0].plot( thetas, data1[:,0], 'b-', label = r"$t$" )
axes[0].plot( thetas, data1[:,1], 'r-' , label = r'$r$')
axes[1].plot( thetas, data2[:,0], 'b' , label = r"$t$" )
axes[1].plot( thetas, data2[:,1], 'r' , label = r"$r$" )

# Draw polarisations
draw_ellipse( 1, 1j , centre = [0.5, 0.9], ax = axes[0], scale=1.2, norm=True, color = 'r' )
draw_ellipse( 1, 3j, centre = [0.5, 0.9], ax = axes[1], scale=1.2, norm=True, color = 'r' )



for theta in thetas[::15]:
    draw_ellipse( np.cos(theta), 1j*np.sin(theta), centre = [theta, -1.1], ax = axes[0], scale=0.12, norm=True, color = 'b' )
    draw_ellipse( np.cos(theta), 1j*np.sin(theta), centre = [theta, -1.1], ax = axes[1], scale=0.12, norm=True, color = 'b' )


axes[0].set_xticks([0, np.pi], ['0', r'$\pi$'])
axes[1].set_xticks([0, np.pi], ['0', r'$\pi$'])
axes[1].set_xlabel( r'$\theta$' )
axes[1].set_ylabel('Amplitude')

axes[0].axis('equal')
axes[1].axis('equal')

axes[0].set_xlim([0, np.pi])
axes[1].set_xlim([0, np.pi])

ax.legend()





#%% Try IXI system (Vogit QD)
S = multilevel_system([10, 10], [0, 0])
S.add_transition( 0, 0, [1, 0] )
S.add_transition( 1, 1, [1, 0] )
S.add_transition( 0, 1, [0, 1j] )
S.add_transition( 1, 0, [0, 1j] )

# At C point
G = WG_Greens_function( [[1, 1j]], [mag], ["WG"])
G.add_mode( [1, 0], 0.1, False, "H loss" )
G.add_mode( [0, 1], 0.1, False, "V loss" )

scat, space = long_photon_scattering(G, S, [1, 0], [1, 0, 0, 0], freq=10)

print('Vogit Geom QD at C point (transmits with phase shift)')
for r, name in zip([complex(i[0,0]) for i in scat], space):
    print( name, " " * (24-len(name)), r)
    
print("Norm = ", np.linalg.norm( np.array(scat) ), "\n" )
# T11 Ended


detunings = np.linspace(0, 2, 100)
data = []

for d in detunings:
    
    S = multilevel_system([10, 10+d], [0, 0])
    S.add_transition( 0, 0, [1, 0] )
    S.add_transition( 1, 1, [1, 0] )
    S.add_transition( 0, 1, [0, 1j] )
    S.add_transition( 1, 0, [0, 1j] )

    # At C point
    G = WG_Greens_function( [[1, 1j]], [mag], ["WG"])
    G.add_mode( [1, 0], 0.1, False, "H loss" )
    G.add_mode( [0, 1], 0.1, False, "V loss" )

    scat, space = long_photon_scattering(G, S, [1, 0], [1, 0, 0, 0], freq=10)
    
    data.append( scat )


#%% Try IXI system (Vogit QD)
# With NON-degenerate ground states
S = multilevel_system([10, 10], [0, 1])
S.add_transition( 0, 0, [1, 0] )
S.add_transition( 1, 1, [1, 0] )
S.add_transition( 0, 1, [0, 1j] )
S.add_transition( 1, 0, [0, 1j] )

# At C point
G = WG_Greens_function( [[1, 1j]], [mag], ["WG"])
G.add_mode( [1, 0], 0.1, False, "H loss" )
G.add_mode( [0, 1], 0.1, False, "V loss" )

scat, space = long_photon_scattering(G, S, [1, 0], [1, 0, 0, 0], freq=10)

print('Vogit Geom QD at C point (transmits with phase shift)')
for r, name in zip([complex(i[0,0]) for i in scat], space):
    print( name, " " * (24-len(name)), r)
    
print("Norm = ", np.linalg.norm( np.array(scat) ), "\n" )
# T12 Ended


#%% Try IXI system (Vogit QD)
# With NON-degenerate ground states
S = multilevel_system([10, 11], [0, 1])
S.add_transition( 0, 0, [1, 0] )
S.add_transition( 1, 1, [1, 0] )
S.add_transition( 0, 1, [0, 1j] )
S.add_transition( 1, 0, [0, 1j] )

# At C point
G = WG_Greens_function( [[1, 1j]], [mag], ["WG"])
G.add_mode( [1, 0], 0.1, False, "H loss" )
G.add_mode( [0, 1], 0.1, False, "V loss" )

scat, space = long_photon_scattering(G, S, [1, 0], [1, 0, 0, 0], freq=10)

print('Vogit Geom QD at C point (transmits with phase shift)')
for r, name in zip([complex(i[0,0]) for i in scat], space):
    print( name, " " * (24-len(name)), r)
    
print("Norm = ", np.linalg.norm( np.array(scat) ), "\n" )
# T12 Ended
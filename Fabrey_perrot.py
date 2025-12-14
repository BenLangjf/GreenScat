# Two atom Fabrey Perrot

# Ben Lang, May 2023

# Idea is to realise a Fabrey Perrot interferometer where the effective spaceing
# between the reflecting atoms is controlled by the rotation of a dipole.


# Referee 2 asked if the transmission matrices of the two individual dipoles could just be multiplied together to get the 2-dipole case. This is an interesting question. One dipole is basically a sort-of-mirror (and/or phase shifter) and the transmission matrix of several mirror/phase shifter’s is just given by multiplying the individiaul transmission matrices. However with two dipoles the excited states hybridise, which is not included by simple multiplication. The symmetric superposition state (superposition of either atom being excited with a “+”sign) and antisymmetric states (“-“ sign) lift apart in energy due to the two dipoles interacting through the waveguide.

# These figures show the complex reflection (colour shade is phase, colour intensity is amplitude), as the distance between the two atoms is changed (x-axis) and the frequency of the input photon is turned. (100 is bang-on resonance). The left figure comes from multiplying the transmission matrices (treating the two dipoles separately), the right figure includes the full ability of the dipoles to hybridise with one another, which introduces some cool looking swirls. Its interesting how the zero reflection points (white lines) kind of all get shifted along one step as they pass the resonance in the full model.

# >>>>> THINK ABOUT Redheffer star product! Might be related.

from MGS_with_position import *


exp = np.exp
cos = np.cos
sin = np.sin
pi = np.pi

# C point waveguide
# This should describe a waveguide with 1 mode in each direction
# at two locations in space. Both C-points.
# Both points experience the same loss into free space, but these losses are uncoupled to one another.

phase_spacing = 0.84757* np.pi  # 1.56745 * np.pi
# this is the exp(ik) phase propegated between the two C-points.

vg = 0.1
mag = (0.5 / vg)

for_mode  = mode([[1,  1j], [1 * np.exp( +1j* phase_spacing ),  1j * np.exp( +1j* phase_spacing )]], mag_list = [np.sqrt(mag), np.sqrt(mag)], direction = +1, name = "WG F")
back_mode = mode([[1, -1j], [1 * np.exp( -1j* phase_spacing ), -1j * np.exp( -1j* phase_spacing )]], mag_list = [np.sqrt(mag), np.sqrt(mag)], direction = -1, name = "WG B")

loss_H_mode_p1 = mode([[1, 0], [1, 0]], mag_list = [np.sqrt(0.1), 0], name = "H loss, p1")
loss_V_mode_p1 = mode([[0, 1], [0, 1]], mag_list = [np.sqrt(0.1), 0], name = "V loss, p1")

loss_H_mode_p2 = mode([[1, 0], [1, 0]], mag_list = [0, np.sqrt(0.1)], name = "H loss, p2")
loss_V_mode_p2 = mode([[0, 1], [0, 1]], mag_list = [0, np.sqrt(0.1)], name = "V loss, p2")

G = WG_Greens_function( [for_mode, back_mode, loss_H_mode_p1, loss_V_mode_p1, loss_H_mode_p2, loss_V_mode_p2 ] )



# This makes the 2 TLS's share a single pair of loss modes.
# Ideally each would have its own set, coupled to the other.


refs = []
trans = []

for theta in np.linspace(0, 2*np.pi, 100):
    # V level system is equivalent to two I systems (in single excitation manifold)
    S = multilevel_system( [100, 100], [0])
    
    S.add_transition( 0, 0, [1, 0], 1, 0 )
    S.add_transition( 1, 0, [cos(theta) , sin(theta) ], 1, 1 ) # At positions 1, not zero!
    
    freq = 100
    
    
    fout, inspace, outspace = S.fanout( G )
    fin  = S.fanin( G )[0]
    
    V = S.full_inner( G )

    
    # Loss modes that couple to each transition SEPERATELY, can be included by adding a little identity to V
    V = V + 0.05 * np.eye(2)
    
    Det = S.detuning_mat( freq ) 
    
    scat_mat =  np.eye(len(inspace)) - fin * lng.inv( V + 1j * Det ) * fout
    instate = np.matrix( [1, 0, 0, 0, 0, 0] ) 
    
    fields = scat_mat * instate.T
    
    refs.append( complex(fields[1]) )
    trans.append( complex(fields[0]) )
    
fig, ax = plt.subplots(figsize = (5,5))
ax.plot( [abs(item) for item in refs], 'r')
ax.plot( [abs(item) for item in trans], 'b')


# Do a series of spectral plots
refs = []
trans = []

fig, ax = plt.subplots(figsize = (5,5))

refs.append( [] )
trans.append( [] )
for freq in np.linspace(90.001, 110, 101):
    
    vg = 0.1
    mag = (0.5 / vg)

    for_mode  = mode([[1,  1j]], mag_list = [np.sqrt(mag), np.sqrt(mag)], direction = +1, name = "WG F")
    back_mode = mode([[1, -1j]], mag_list = [np.sqrt(mag), np.sqrt(mag)], direction = -1, name = "WG B")

    loss_H_mode_p1 = mode([[1, 0], [1, 0]], mag_list = [np.sqrt(0.1), 0], name = "H loss, p1")
    loss_V_mode_p1 = mode([[0, 1], [0, 1]], mag_list = [np.sqrt(0.1), 0], name = "V loss, p1")

    G = WG_Greens_function( [for_mode, back_mode, loss_H_mode_p1, loss_V_mode_p1 ] )

    S = multilevel_system( [100], [0])
    S.add_transition( 0, 0, [1, 0], 1, 0 )
    
    fout, inspace, outspace = S.fanout( G )
    fin  = S.fanin( G )[0]
    
    V = S.full_inner( G )
    
    Det = S.detuning_mat( freq ) 
    
    scat_mat =  np.eye(len(inspace)) - fin * lng.inv( V + 1j * Det ) * fout
    instate = np.matrix( [1, 0, 0, 0] ) 
    
    fields = scat_mat * instate.T
    
    refs[-1].append( complex(fields[1]) )
    trans[-1].append( complex(fields[0]) )


ax.plot( np.linspace(90.001, 110, 101)  - 100,  [abs(item) for item in refs[-1]], 'g')


def add_new_plot( theta, distance, ax, line = "r" ):
    
    refs = []
    trans = []
    Vs = []
    for freq in np.linspace(90.1, 110, 401):
        
        phase_spacing = distance * 2 * pi * (freq/100)
        
        vg = 0.1
        mag = (0.5 / vg)
    
        for_mode  = mode([[1,  1j], [1 * np.exp( +1j* phase_spacing ),  1j * np.exp( +1j* phase_spacing )]], mag_list = [np.sqrt(mag), np.sqrt(mag)], direction = +1, name = "WG F")
        back_mode = mode([[1, -1j], [1 * np.exp( -1j* phase_spacing ), -1j * np.exp( -1j* phase_spacing )]], mag_list = [np.sqrt(mag), np.sqrt(mag)], direction = -1, name = "WG B")
    
        loss_H_mode_p1 = mode([[1, 0], [1, 0]], mag_list = [np.sqrt(0.1), 0], name = "H loss, p1")
        loss_V_mode_p1 = mode([[0, 1], [0, 1]], mag_list = [np.sqrt(0.1), 0], name = "V loss, p1")
    
        loss_H_mode_p2 = mode([[1, 0], [1, 0]], mag_list = [0, np.sqrt(0.1)], name = "H loss, p2")
        loss_V_mode_p2 = mode([[0, 1], [0, 1]], mag_list = [0, np.sqrt(0.1)], name = "V loss, p2")
    
        G = WG_Greens_function( [for_mode, back_mode, loss_H_mode_p1, loss_V_mode_p1, loss_H_mode_p2, loss_V_mode_p2 ] )
    
        
        # V level system is equivalent to two I systems (in single excitation manifold)
        S = multilevel_system( [100, 100], [0])
        
        S.add_transition( 0, 0, [1, 0], 1, 0 )
        S.add_transition( 1, 0, [cos(theta) , sin(theta) ], 1, 1 ) # At positions 1, not zero!
        
        fout, inspace, outspace = S.fanout( G )
        fin  = S.fanin( G )[0]
        
        V = S.full_inner( G )
        
        Det = S.detuning_mat( freq ) 
        
        scat_mat =  np.eye(len(inspace)) - fin * lng.inv( V + 1j * Det ) * fout
        instate = np.matrix( [1, 0, 0, 0, 0, 0] ) 
        
        fields = scat_mat * instate.T
        
        refs.append( complex(fields[1]) )
        trans.append( complex(fields[0]) )
        Vs.append( V )
    
    
    ax.plot( np.linspace(90.001, 110, 401) - 100, [abs(item) for item in refs], line)
    ax.set_xlabel(r'Detuning ($\delta$)')
    ax.set_ylabel(r'$|r|$')
    
    return refs, trans, Vs


add_new_plot(0, 9, ax, 'r')
add_new_plot(0, 9.25, ax, 'b')
add_new_plot(0.5*pi, 9.25, ax, 'k')

# for theta, distance in [[0, 18], [0, 18.64], [0.3*pi, 18.64], [0.1*pi, 18.64], [0.5*pi, 18] ]: #, [0, 2.25], [0, 2.5], [0.5*pi, 2.5]]:
    
#     # Do a series of spectral plots
#     refs.append( [] )
#     trans.append( [] )
    
#     for freq in np.linspace(90.001, 110, 401):
        
#         phase_spacing = distance * pi * (100/freq)
        
#         vg = 0.1
#         mag = (0.5 / vg)
    
#         for_mode  = mode([[1,  1j], [1 * np.exp( +1j* phase_spacing ),  1j * np.exp( +1j* phase_spacing )]], mag_list = [np.sqrt(mag), np.sqrt(mag)], name = "WG F")
#         back_mode = mode([[1, -1j], [1 * np.exp( -1j* phase_spacing ), -1j * np.exp( -1j* phase_spacing )]], mag_list = [np.sqrt(mag), np.sqrt(mag)], name = "WG B")
    
#         loss_H_mode_p1 = mode([[1, 0], [1, 0]], mag_list = [np.sqrt(0.1), 0], name = "H loss, p1")
#         loss_V_mode_p1 = mode([[0, 1], [0, 1]], mag_list = [np.sqrt(0.1), 0], name = "V loss, p1")
    
#         loss_H_mode_p2 = mode([[1, 0], [1, 0]], mag_list = [0, np.sqrt(0.1)], name = "H loss, p2")
#         loss_V_mode_p2 = mode([[0, 1], [0, 1]], mag_list = [0, np.sqrt(0.1)], name = "V loss, p2")
    
#         G = WG_Greens_function( [for_mode, back_mode, loss_H_mode_p1, loss_V_mode_p1, loss_H_mode_p2, loss_V_mode_p2 ] )
    
        
#         # V level system is equivalent to two I systems (in single excitation manifold)
#         S = multilevel_system( [100, 100], [0])
        
#         S.add_transition( 0, 0, [1, 0], 1, 0 )
#         S.add_transition( 1, 0, [cos(theta) , sin(theta) ], 1, 1 ) # At positions 1, not zero!
        
#         fout, inspace, outspace = S.fanout( G )
#         fin  = S.fanin( G )[0]
        
#         V = S.full_inner( G )
        
#         Det = S.detuning_mat( freq ) 
        
#         scat_mat =  np.eye(len(inspace)) - fin * lng.inv( V + 1j * Det ) * fout
#         instate = np.matrix( [1, 0, 0, 0, 0, 0] ) 
        
#         fields = scat_mat * instate.T
        
#         refs[-1].append( complex(fields[1]) )
#         trans[-1].append( complex(fields[0]) )
    

#     ax.plot( np.linspace(90.001, 110, 401) - 100, [abs(item) for item in refs[-1]], 'r')

for distance, theta in [[18, 0], [18.5,0], [18.5,0.5*pi], [18.3,0], [18.4,0]]:
    fig, ax = plt.subplots(figsize = (5,5))
    add_new_plot(theta, distance, ax)
    ax.set_xlabel(r'Detuning ($\delta$)')
    ax.set_ylabel(r'$|r|$')



    



### Colour (2D) Plot

from phase_amp import *

fig, axes = plt.subplots(5, 1, figsize = (5,16))
n=0

for theta in [0, 0.125*pi, 0.25*pi, 0.375*pi, 0.5*pi]:
    
    refs = []
    trans = []
    Vs = []
    
    distances = np.linspace(0, 7.5, 300)
    freqs = np.linspace(90.001, 110, 201)
    
    for distance in distances:
        
        # Do a series of spectral plots
        refs.append( [] )
        trans.append( [] )
        Vs.append( [] )
        
        for freq in freqs:
            
            phase_spacing = distance * 2 * pi * (freq/100)
            
            vg = 0.2 #0.1
            mag = (0.5 / vg)
            loss_strength = 0.3 #0.1
            
            half_phase = np.exp( +1j * phase_spacing * 0.5 )
            conj = np.conj
            
            for_mode  = mode([[1 * conj(half_phase),  1j * conj(half_phase)], [1 * half_phase,  1j * half_phase]], mag_list = [np.sqrt(mag), np.sqrt(mag)], direction = +1, name = "WG F")
            back_mode = mode([[1 * half_phase, -1j * half_phase], [1 * conj(half_phase), -1j * conj(half_phase)]], mag_list = [np.sqrt(mag), np.sqrt(mag)], direction = -1, name = "WG B")
        
            loss_H_mode_p1 = mode([[1, 0], [1, 0]], mag_list = [np.sqrt(loss_strength), 0], name = "H loss, p1")
            loss_V_mode_p1 = mode([[0, 1], [0, 1]], mag_list = [np.sqrt(loss_strength), 0], name = "V loss, p1")
        
            loss_H_mode_p2 = mode([[1, 0], [1, 0]], mag_list = [0, np.sqrt(loss_strength)], name = "H loss, p2")
            loss_V_mode_p2 = mode([[0, 1], [0, 1]], mag_list = [0, np.sqrt(loss_strength)], name = "V loss, p2")
        
            G = WG_Greens_function( [for_mode, back_mode, loss_H_mode_p1, loss_V_mode_p1, loss_H_mode_p2, loss_V_mode_p2 ] )
        
            
            # V level system is equivalent to two I systems (in single excitation manifold)
            S = multilevel_system( [100, 100], [0])
            
            S.add_transition( 0, 0, [1, 0], 1, 0 )
            S.add_transition( 1, 0, [cos(theta) , sin(theta) ], 1, 1 ) # At positions 1, not zero!
            
            fout, inspace, outspace = S.fanout( G )
            fin  = S.fanin( G )[0]
            
            V = S.full_inner( G )
            
            Det = S.detuning_mat( freq ) 
            
            scat_mat =  np.eye(len(inspace)) - fin * lng.inv( V + 1j * Det ) * fout
            instate = np.matrix( [1, 0, 0, 0, 0, 0] ) 
            
            fields = scat_mat * instate.T
            
            refs[-1].append( complex(fields[1]) )
            trans[-1].append( complex(fields[0]) )
            Vs[-1].append( V )

    refs = np.array( refs )
    trans = np.array( trans )
    
    phase_amp( refs.transpose(), distances, np.linspace(90.001, 110, 201), axes[n], aspect = 0.25 )
    n +=1 



### Important Question
# Is my move complicated "V" system model equivalent to just multiplying the transmision matrices?

# 1 atom t-matrix as a function of freqs
mats = []
freqs = np.linspace(90.001, 110, 201)
 
for freq in freqs:
    
    # vg = 0.1
    # mag = (0.5 / vg)

    for_mode  = mode([[1,  1j]], mag_list = [np.sqrt(mag)], direction = +1, name = "WG F")
    back_mode = mode([[1, -1j]], mag_list = [np.sqrt(mag)], direction = -1, name = "WG B")

    loss_H_mode_p1 = mode([[1, 0], [1, 0]], mag_list = [np.sqrt(loss_strength), 0], name = "H loss, p1")
    loss_V_mode_p1 = mode([[0, 1], [0, 1]], mag_list = [np.sqrt(loss_strength), 0], name = "V loss, p1")

    G = WG_Greens_function( [for_mode, back_mode, loss_H_mode_p1, loss_V_mode_p1 ] )

    S = multilevel_system( [100], [0])
    S.add_transition( 0, 0, [1, 0], 1, 0 )
    
    fout, inspace, outspace = S.fanout( G )
    fin  = S.fanin( G )[0]
    
    V = S.full_inner( G )
    
    Det = S.detuning_mat( freq ) 
    
    scat_mat =  np.eye(len(inspace)) - fin * lng.inv( V + 1j * Det ) * fout
    mats.append( scat_mat )


distances = np.linspace(0, 7.5, 300)

# Now to check multiplicity
refs_mul = []

for distance in distances:
    
    # Do a series of spectral plots
    refs_mul.append( [] )
    
    for one_atom_mat, freq in zip(mats, freqs):
        
        k = 2* pi * (freq/100)
        
        propegation_mat = np.matrix( [[ np.exp(1j*distance*k), 0, 0, 0 ], [ 0, np.exp(-1j*distance*k), 0, 0 ], [0, 0, 0, 0], [0, 0, 0, 0]] )
        
        tot_mat = one_atom_mat * propegation_mat * one_atom_mat
        
        instate = np.matrix( [1, 0, 0, 0] )
        fields = tot_mat * instate.T
        
        refs_mul[-1].append( complex(fields[1]) )
    
refs_mul = np.array( refs_mul )
phase_amp( refs_mul.transpose(), distances, np.linspace(90.001, 110, 201), aspect = 0.25 )
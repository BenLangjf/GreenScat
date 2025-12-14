
from MGE import *


#%% Emission Tests

# Keep Group velocity fixed throughout
vg = 0.1
mag = (0.5 / vg)

# Test 1, decay of a two-level-system
G = WG_Greens_function( [], [], [],)
G.add_mode( [1,1], mag, include_back_partner = False )

S = multilevel_system([1000], [0])
S.add_transition( 0, 0, [1,  0] )

_, labs_g, labs_e = S.fanout( G )
full_labs = labs_e + labs_g

H_size = len(full_labs)


L, H = decay_from_excited_state( G, S )
res = mesolve(L, basis(H_size, 0), np.linspace(0, 3, 300))

c = 'rgkbyc' * 5

fig, ax = plt.subplots(figsize = (5,5))

for n in range(H_size):
    op_list = [0]*H_size
    op_list[n] = 1
    
    ax.plot( np.linspace(0, 3, 300), [(op * Qobj(np.diag(op_list)) ).tr() for op in res.states]
            , c[n] + 'x' )
ax.legend()
ax.set_title('Basic Decay')
ax.set_xlim([0, 3])
ax.set_ylim([0, 1])

#%% Test 2.1, matching Stephen Hughes
# Fig 2(a) in https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.118.063601
G = WG_Greens_function( [[1, 1]], [mag], ["WG"])

S = multilevel_system([1000, 1000], [0])
S.add_transition( 0, 0, [1,  1j] )
S.add_transition( 1, 0, [1, -1j] )

_, labs_g, labs_e = S.fanout( G )
full_labs = labs_e + labs_g

H_size = len(full_labs)

L, H = decay_from_excited_state( G, S )
res = mesolve(L, basis(H_size, 0), np.linspace(0, 3, 300))

c = 'rgkbyc' * 5

fig, ax = plt.subplots(figsize = (5,5))
# This *should* match figure 2(a) in https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.118.063601
# However, it does NOT.
# I beleive the linked paper has ommited a factor of "i". The delta term in the linked paper's equation (3)
# has a factor of i that I beleive they accidentally neglected from their calculations. If I leave
# out this factor of i, then I reproduce the graph perfectly.
# This factor of i has no effect at all on most of the paper's figures. It makes 2(a) look
# slightly different, but not in a way that changes the paper's conclusions.


for n in range(H_size):
    
    ax.plot( np.linspace(0, 3, 300), [complex(op[n,n]) for op in res.states]
            , c[n] + 'x', label=full_labs[n] )
ax.legend()
ax.set_title('Matching Stephen Hughes (2a)')
ax.set_xlim([0, 3])
ax.set_ylim([0, 1])


# Fig 2(b)
G = WG_Greens_function( [[1, 0]], [mag], ["WG"])

S = multilevel_system([1000, 1000], [0])
S.add_transition( 0, 0, [1,  1j] )
S.add_transition( 1, 0, [1, -1j] )

_, labs_g, labs_e = S.fanout( G )
full_labs = labs_e + labs_g

H_size = len(full_labs)

L, H = decay_from_excited_state( G, S )
res = mesolve(L, basis(H_size, 0), np.linspace(0, 3, 300))

c = 'rgkbyc' * 5

fig, ax = plt.subplots(figsize = (5,5))
# This should match figure 2(b) in https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.118.063601

for n in range(H_size):
    
    ax.plot( np.linspace(0, 3, 300), [complex(op[n,n]) for op in res.states]
            , c[n] + 'x', label=full_labs[n] )
ax.legend()
ax.set_title('Matching Stephen Hughes (2b)')
ax.set_xlim([0, 3])
ax.set_ylim([0, 1])

        
#%% Directionality, directly
G = WG_Greens_function( [[1, 1j]], [mag], ["WG"])

S = multilevel_system([1000, 1000], [0])
S.add_transition( 0, 0, [1,  1j] )
S.add_transition( 1, 0, [1, -1j] )

_, labs_g, labs_e = S.fanout( G )
full_labs = labs_e + labs_g

H_size = len(full_labs)

L, H = decay_from_excited_state( G, S )
res = mesolve(L, basis(H_size, 0), np.linspace(0, 3, 300))

c = 'rgkbyc' * 5

fig, ax = plt.subplots(figsize = (5,5))
# This should match figure 2(a) in https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.118.063601

for n in range(H_size):
    
    ax.plot( np.linspace(0, 3, 300), [complex(op[n,n]) for op in res.states]
            , c[n] + 'x', label=full_labs[n] )
ax.legend()
ax.set_title('Normal Directionality')
ax.set_xlim([0, 3])
ax.set_ylim([0, 1])

#%% Directionality, via superposition
G = WG_Greens_function( [[1, 1j]], [mag], ["WG"])

S = multilevel_system([1000, 1000], [0])
S.add_transition( 0, 0, [1, 0] )
S.add_transition( 1, 0, [0, 1] )

_, labs_g, labs_e = S.fanout( G )
full_labs = labs_e + labs_g

H_size = len(full_labs)

L, H = decay_from_excited_state( G, S )

# Consider superposition excited state (aiming for directionality)
res = mesolve(L, (basis(H_size, 0) - 1j * basis(H_size, 1))/np.sqrt(2), np.linspace(0, 3, 300))

c = 'rgkbyc' * 5

fig, ax = plt.subplots(figsize = (5,5))

for n in range(H_size):
    op_list = [0]*H_size
    op_list[n] = 1
    
    ax.plot( np.linspace(0, 3, 300), [(op * Qobj(np.diag(op_list)) ).tr() for op in res.states]
            , c[n] + 'x', label=full_labs[n] )
ax.legend()
ax.set_title('Superposition Directionality')
ax.set_xlim([0, 3])
ax.set_ylim([0, 1])

#%% Directionality, via superposition at an elliptical point
# This situation seems paradoxical, directionality shouldn't work because
# that would allow non-orthogonal quantum states to be resolved perfectly.
# But..., the simplification of the model would imply this does happen.
# So what extra details of the proper model prevent it?
G = WG_Greens_function( [[2, 1j]], [mag], ["WG"])

S = multilevel_system([1000, 1000], [0])
S.add_transition( 0, 0, [1, 0] )
S.add_transition( 1, 0, [0, 1] )

_, labs_g, labs_e = S.fanout( G )
full_labs = labs_e + labs_g

H_size = len(full_labs)

L, H = decay_from_excited_state( G, S )

# Consider superposition excited state (aiming for directionality)
res1 = mesolve(L, (basis(H_size, 0) - 2j * basis(H_size, 1))/np.sqrt(5), np.linspace(0, 5, 500))

c = 'rgkbyc' * 5

fig, ax = plt.subplots(figsize = (5,5))

for n in range(H_size):
    op_list = [0]*H_size
    op_list[n] = 1
    
    ax.plot( np.linspace(0, 5, 500), [(op * Qobj(np.diag(op_list)) ).tr() for op in res1.states]
            , c[n] + '-', label=full_labs[n] )
ax.legend()
ax.set_title('Ellipse Directionality Paradox')
ax.set_xlim([0, 3])
ax.set_ylim([0, 1])

fig.savefig("Ellipe_lambda_paradox_resolution.pdf", bbox_inches='tight')

ax = S.represent()
draw_ellipse( 2, 1j, centre = [10, 10], ax = ax, scale=28, norm=True )

fig = ax.get_figure()
fig.savefig("lambda_system.pdf", bbox_inches='tight')


#%% Directionality, via superposition at an elliptical point
# This time I will have a state that isn't fully directional immediately, but does that at later times.
G = WG_Greens_function( [[2, 1j]], [mag], ["WG"])

S = multilevel_system([1000, 1000], [0])
S.add_transition( 0, 0, [1, 0] )
S.add_transition( 1, 0, [0, 1] )

_, labs_g, labs_e = S.fanout( G )
full_labs = labs_e + labs_g

H_size = len(full_labs)

L, H = decay_from_excited_state( G, S )

# Consider superposition excited state (aiming for directionality)
res3 = mesolve(L, (basis(H_size, 0) - 1j * basis(H_size, 1))/np.sqrt(2), np.linspace(0, 5, 500))

c = 'rgkbyc' * 5

fig, ax = plt.subplots(figsize = (5,5))

for n in range(H_size):
    op_list = [0]*H_size
    op_list[n] = 1
    
    ax.plot( np.linspace(0, 5, 500), [(op * Qobj(np.diag(op_list)) ).tr() for op in res3.states]
            , c[n] + '-', label=full_labs[n] )
ax.legend()
ax.set_title('Ellipse Directionality Paradox, later times')
ax.set_xlim([0, 3])
ax.set_ylim([0, 1])

fig.savefig("Ellipe_lambda_paradox_resolution_directional_later.pdf", bbox_inches='tight')



#%% This situation seems paradoxical, directionality shouldn't work because
# that would allow non-orthogonal quantum states to be resolved perfectly.
# But..., the simplification of the model would imply this does happen.
# So what extra details of the proper model prevent it?
G = WG_Greens_function( [[2, 1j]], [mag], ["WG"])

S = multilevel_system([1000, 1000], [0])
S.add_transition( 0, 0, [1, 0] )
S.add_transition( 1, 0, [0, 1] )

_, labs_g, labs_e = S.fanout( G )
full_labs = labs_e + labs_g

H_size = len(full_labs)

L, H = decay_from_excited_state( G, S )

# Consider superposition excited state (aiming for directionality)
res2 = mesolve(L, (basis(H_size, 0) + 2j * basis(H_size, 1))/np.sqrt(5), np.linspace(0, 5, 500))

c = 'rgkbyc' * 5

fig, ax = plt.subplots(figsize = (5,5))

for n in range(H_size):
    op_list = [0]*H_size
    op_list[n] = 1
    
    ax.plot( np.linspace(0, 5, 500), [(op * Qobj(np.diag(op_list)) ).tr() for op in res2.states]
            , c[n] + 'x', label=full_labs[n] )
ax.legend()
ax.set_title('Ellipse Directionality Paradox Conjuate')
ax.set_xlim([0, 3])
ax.set_ylim([0, 1])

# if I consider time another degree of freedom, then its like I have a bigger density matrix to account for this.
# This is done in the most obvious way, slicing up the emission process into a series of times-steps and having the time
# frame of the emission act like a degree of freedom.
previous_photon_density_matrix_1 = Qobj(np.zeros((2,2)))
previous_photon_density_matrix_2 = Qobj(np.zeros((2,2)))

matrix_stream_1 = []
matrix_stream_2 = []

for s1, s2 in zip(res1.states[1:], res2.states[1:]):
    # Post-selecting requires me to scale this up...
    # Get a photon-only density matrix
    photon_density_matrix_1 = Qobj(s1.full()[2:,2:])
    photon_density_matrix_2 = Qobj(s2.full()[2:,2:])
    
    marginal_density_1 = photon_density_matrix_1 - previous_photon_density_matrix_1
    marginal_density_2 = photon_density_matrix_2 - previous_photon_density_matrix_2
    
    previous_photon_density_matrix_1 = photon_density_matrix_1
    previous_photon_density_matrix_2 = photon_density_matrix_2
    
    matrix_stream_1.append( marginal_density_1 )
    matrix_stream_2.append( marginal_density_2 )

# I don't think I can just take an incoherent sum of all the matrices in each steam (IE put them diagonal)
# But lets take this as a first simple thing.

L = len(matrix_stream_1)
zero = 0 * photon_density_matrix_1.full()

# Assume coherences between the waves at different times.
# These must exist, as the light wave has a wavelength.
# Note, that the phases of these coherences would in real life exist, but are not relevant here.
total_density_over_time_1 = Qobj( np.block( [ [ np.sqrt(item1.full())*np.sqrt(item2.full()) for item1 in matrix_stream_1]
                                             for item2 in matrix_stream_1 ] ) )

total_density_over_time_2 = Qobj( np.block( [ [ np.sqrt(item1.full())*np.sqrt(item2.full()) for item1 in matrix_stream_2]
                                             for item2 in matrix_stream_2 ] ) )

print('Purity, should match purity of final state',
      (photon_density_matrix_1**2).tr(),
      (total_density_over_time_1**2).tr(),
      (total_density_over_time_2**2).tr() )


fig, ax = plt.subplots( figsize = (5,5))
ax.plot( np.diag( total_density_over_time_1.full())[::2] )
ax.plot( np.diag( total_density_over_time_1.full())[1::2] )

print('trace 1', total_density_over_time_1.tr() )
print('trace 2', total_density_over_time_2.tr() )

print('Overlap from Fiedlity', fidelity(total_density_over_time_1,total_density_over_time_2)**2 )
print('Overlap from Tr_dist', (1 - tracedist(total_density_over_time_1,total_density_over_time_2)**2) )
print('Expected overlap', fidelity(res1.states[0], res2.states[0])**2, (1 - tracedist(res1.states[0], res2.states[0])**2) )

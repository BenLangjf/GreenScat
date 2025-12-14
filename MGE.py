# MGS: Multilevel system Greens Emission.

# Calculate emission of multi-level quantum systems coupled to optical enviroments

# Ben Lang 2023

# Notes:
# I have Hilbert space, hosting quamtum wavefunctions. The operators in this space are matrices.
# I also have real space, containing complex vectors with an inner product  a dot b*
# The entires in the operators are determined by these real space inner products.

# Non-diagonal disipators: https://arxiv.org/pdf/2111.04041
# https://en.wikipedia.org/wiki/Redheffer_star_product


from MGS import *
from qutip import *


def decay_from_excited_state( G, S ):
    """
    

    Parameters
    ----------
    G : Greens function object
        The Green's function describing the optical modes.
    S : Multilevel system object
        The Multilevel quantum system.

    Returns
    -------
    Liouvillian : QuTiP Superoperator
        The Liouvillian decribing time evoltion of the open system.
    Hamiltonian : QuTiP Operator
        The system Hamiltonian.

    """
    # Decay from excited state, based on my model.
    # find sizes of excited and ground state manifolds.
    
    # Freqs must include energies of both excited states and ground states;
    # Freqs = [E_1, E_2,..., G_1, G_2 ....]
    
    e_man = S.N_ext
    g_man = S.N_grd * len(G.G_list)
    
    
    # Note g_man has size: number of modes times number of ground states.
    full_man = e_man + g_man

    Hamil = 0
    Liouv = 0
    
    hbar = 1
    
    def dec_op(x, n):
        # First e_man states are excited states
        # rest are ground tensor photon states.
        # Form |g_n><e_x|
        return basis(full_man, e_man + n) * basis(full_man, x).dag()
    

    # Frequency (energy) Hamiltonian terms (diagonals)
    # kronecker sum of light and ground state energies
    energy = S.ext_ens + S.grd_ens * len(G.G_list)
    
    
    Hamil += Qobj( np.diag(energy) )
    
    Liouv += liouvillian(Hamil)  #(-1j / hbar) *(spre(Hamil) - spost(Hamil)) 
        
    # W depends on two choices of excited state, and two of ground state.
    W = np.zeros( (e_man, e_man, g_man, g_man ), dtype = complex )
    
    # Coupling/disipation terms
    for n in range(g_man):
        for m in range(g_man):
            
            # n indicates a combination of optical mode and ground state.
            g_n, mode_n = divmod( n, len(G.name_list) )
            g_m, mode_m = divmod( m, len(G.name_list) )
                    
            W[:, :, n, m] =  1j * S.partial_inner( G, g_n, g_m, mode_n, mode_m )
                    
    # Gamma is similar to W, but it summs over the ground states, and only the diagonal ones
    Gamma = 0 + 0j
    for n in range(g_man):
        Gamma += W[:, :, n, n]
       
    ## Now settup Liov and Hamil
    
    # I am using the "General Form" where the operators sandwiching the state can be different.
    # Diagonalisation can be used to recude to a speial form where you can make the 
    # operators either side symmetric. But this is not being done here.
    
    # Sandwich terms
    for x in range(e_man):
        for y in range(e_man):
            for n in range(g_man):
                for m in range(g_man):
                    # |g_n><e_x| rho |e_y><g_m|
                    a = dec_op(x, n)
                    b = dec_op(y, m)

                    # Note, the x-y and n-m indices swap for the second W.                    
                    Liouv += spre(a) * spost(b.dag()) * ( W[x, y, n, m] - np.conj( W[y, x, m, n] ) ) * (-1j)
    
    ## Non-Sandwich terms:
    # New Derivatioon
    # This seems to either be very almost right, or actually right
    # It is also in agreement with the equation in Stephen Hughes
    # (although I include real diagonals which he excldes)
    # and I am mis-reading the figures in comparisons to Stephen Hughes.
    for x in range(e_man):
        for y in range(e_man):
            # |e_x><e_y| rho
            # and rho |e_y><e_x|
            Opp = basis(full_man, x) * basis(full_man, y).dag()
              
            # Note that xy swap for the second Gamma.
            Liouv += spre(Opp)  * +1j * Gamma[x, y]
            Liouv += spost(Opp) * -1j * np.conj(Gamma[y, x])
    
    # # Non-Sandwich terms (old):
    # To Match the paper of Stephen Hughes (https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.118.063601)
    # we can instead use thee non-sandwich terms. These basically drop a factor of "i" from the real part terms.
    # I beleive this is a mistake that the paper made so this version matches its graphs.
    # However, the version above is correct, so should be used.
    # for x in range(e_man):
    #     for y in range(e_man):
    #         # |e_x><e_y| rho
    #         # and rho |e_y><e_x|
    #         Opp = basis(full_man, x) * basis(full_man, y).dag()
                       
    #         Liouv += ( spre(Opp) + spost(Opp) ) * np.imag(Gamma[x, y]) * (-1)  # hbar = 1 untis.
    #         Liouv += ( spre(Opp) - spost(Opp) ) * np.real(Gamma[x, y]) * (1) # Do NOT multiply by 1j.

    return Liouv, Hamil   #, Hamil # Returns the Liovillian and Hamiltonian

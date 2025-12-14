# MGS: Multilevel system Greens Scattering tool.

# Calculate scattering and emission of multi-level quantum systems coupled to optical enviroments

# Ben Lang 2023

# Notes:
# I have Hilbert space, hosting quamtum wavefunctions. The operators in this space are matrices.
# I also have real space, containing complex vectors with an inner product  a dot b*
# The entires in the operators are determined by these real space inner products.

# Non-diagonal disipators: https://arxiv.org/pdf/2111.04041
# https://en.wikipedia.org/wiki/Redheffer_star_product


import numpy as np
import numpy.linalg as lng
import matplotlib.pyplot as plt
from scipy.linalg import eig as eigs


def space_kron( names_1, names_2 ):
    """
    Creates a new name list, giving all combinations of the two lists provided.

    Parameters
    ----------
    names_1 : list of strings
        A list of mode names.
    names_1 : list of strings
        A list of mode names.
    
    Returns
    -------
    total_space : list of strings
    """
    
    total_space = []
    
    for item1 in names_1:
        for item2 in names_2:
            total_space.append( item1 + ", " + item2 )
            
    
    return total_space



def draw_ellipse( Ex, Ey, centre = [0, 0], ax = "none", scale=0.8, norm=True, color = 'k' ):
    """
    Draw the polarisation ellipse defined by Ex, Ey

    Parameters
    ----------
    Ex : complex
        The X compoenent of the electric field.
    Ey : complex
        The Y compoenent of the electric field.
    centre : list
        centre point coordinates of the ellipse.
    ax : matplotlib axis or "none"
        the axis to draw the ellipse on. For "none" makes a new one.
    scale : float
        scale of ellipse
    norm : bool
        Whether to normalise the given Ex, Ey first, or keep its scale.
    color : string
        colour of the ellipse

    """
    
    if ax == "none":
        fig, ax = plt.subplots(figsize = (5,5))
    
    t_range = np.linspace(0, 2*np.pi, 40)
    
    # Normalise:
    if norm:
        N = 2*np.linalg.norm([Ex, Ey])
        nEx = scale * Ex / N
        nEy = scale * Ey / N
    else:
        nEx = Ex
        nEy = Ey
    
    Xvals = np.real( nEx * np.exp(1j*t_range) )
    Yvals = np.real( nEy * np.exp(1j*t_range) )
    
    ## Add arrowhead:
        
    ax.plot( [centre[0] + x for x in Xvals], [centre[1] + y for y in Yvals], color + '-' )
    
    ## Add arrowhead:
    ax.arrow( centre[0] + Xvals[-2], centre[1] + Yvals[-2], Xvals[-1] - Xvals[-2], Yvals[-1] - Yvals[-2], 
             shape='full', lw=0, length_includes_head=True, head_width=.25*scale, color = color)

class WG_Greens_function:
    "Waveguide Greens Function Class"
    # This class is a representation of the greens function of a waveguide.
    # Multimode waveguides _should_ be possibly by providing a list of polarisations and group velocities (not tested)
    
    def __init__(self, polarisation_list = [[1, 1j]], mag_list = [1], names = [] ):
        # Initialised using 2D vector for the polarisation of the WG at the atom location.
        # and the group velocity (in units of c).
        # zeta (for a WG mode) is: zeta = a / (2*vg)
        # For non WG modes put other zetas.
        
        self.G_list = []
        self.name_list = []
        self.E_list = []
        
        while len(names)< len(polarisation_list):
            names.append('No name')
        
        for pol, mag, name in zip( polarisation_list, mag_list, names ):
            self.add_mode( pol, mag, True, name )
            
    
    def add_mode(self, polarisation, mag=1.0, include_back_partner = True, name = "No name" ):
        """
        Add a new electromagnetic mode to the Greens function (or pair of modes).

        Parameters
        ----------
        polarisation : 2 element list of complex
            The polarisation of the new mode, at the location the Greens function describes.
        mag : float
            The "strength" (density of states) of the mode at the location.
        include_back_partner : bool
            If true, both a forward and backward mode are added, backward having
            the conjugate polarisation.
        name : str
            name of the mode.
        """
        
        # Normalise polarisation vector to length 1.
        polarisation = np.sqrt(mag) * np.array(polarisation) / np.linalg.norm( polarisation )
        
        E_f = np.matrix( polarisation )
        
        self.G_list.append( 0.5 * E_f.H * E_f  )
        self.E_list.append( E_f )
       
        if include_back_partner:
            self.name_list.append( name + "-for" )
            
            E_b = np.matrix( np.conj(polarisation) )
            self.G_list.append( 0.5 * E_b.H * E_b )
            self.E_list.append( E_b )
        
            self.name_list.append( name + "-back" )
        else:
            self.name_list.append( name )
    
    
class multilevel_system:
    "Multilevel System Class"
    
    def __init__( self, ext_ens, grd_ens, e_names=False, g_names=False ):
        # System has some number of excited and ground states.
        # dipoles should be an N_ext by N_grd array, with each element a 2-vector.
        
        N_ext, N_grd = len(ext_ens), len(grd_ens)
        
        self.N_ext = N_ext
        self.N_grd = N_grd
        
        if e_names:
            self.e_names = e_names
        else:
            self.e_names = ["e" + str(n+1) for n in range(N_ext)]
        
        
        if g_names:
            self.g_names = g_names
        else:
            self.g_names = ["g" + str(n+1) for n in range(N_grd)]
            
        # Initialise to zero dipoles.
        self.dipoles = [ [ np.matrix([0, 0]) for _ in range(N_grd) ] for _ in range(N_ext) ]
                
        # Set Energy levels
        self.ext_ens = ext_ens
        self.grd_ens = grd_ens
        
        
    def add_transition( self, en, gn, dipole, normalise = True, omega_scale = False ):
        """
        Add a new dipole allowed transition.

        Parameters
        ----------
        en : int
            The excited state involved
        gn : int
            The ground state involved
        dipole : 2 element list of complex
            The dipole of the transition.
        normalise : bool
            If true, the dipole is normalised to unit length.
        omega_scale : bool
            If true, the dipole magnitude is scaled by sqrt(omega)
        """
        
        # Normalise dipole vector to length 1.
        if normalise:
            dipole = dipole / np.linalg.norm( dipole )
            
        # In a waveguide, coupling strength is proportional to sqrt(omega)
        # (In other types of optical mode its proportional to other functions of omega)
        # If desired, this can be accoutned for by scaling the dipole scale by sqrt(omega)
        # Note that, if this is used in a scattering problem where the input photon
        # is significantly deturned from the transition then its not exactly right, 
        # as the photon frequency should be used. However, for strongly detuned photons,
        # the interaction strength is near zero anyway, so ultimately I don't 
        # beleive the approximation posses any issue.
        
        if omega_scale:
            dipole = dipole * np.sqrt( self.ext_ens[en] - self.grd_ens[gn] )
        
        self.dipoles[en][gn] = np.matrix( dipole )
        
    
    def fanout( self, Greens ):
        """
        Procudes a matrix that maps from states where the excitation is photonic
        to those where it is material.
        So maps from an input matrix of amplitudes of form ['WG-for, g1','WG-for, g2', 'WG-back, g1' ...
        To one like [ e1, e2, e3, ... ]
        """
        
        gsize = len(self.g_names)
        esize = len(self.e_names)
        Esize = len(Greens.name_list)
        
        in_space  = space_kron( Greens.name_list, self.g_names )
        out_space = self.e_names
        
        out_mat = np.zeros( ( len(out_space), len(in_space) ), dtype= complex )
        
        for j1 in range( Esize ):
            for j2 in range( gsize ):
                for i in range( esize ):
                    j = gsize*j1 + j2
                    
                    # .item() converts 1x1 matrix into complex scalar.
                    out_mat[i, j] = ( Greens.E_list[ j1 ] * self.dipoles[i][j2].H ).item() # d* dot E
                    
        return np.matrix( out_mat ), in_space, out_space
    
    
    def fanin( self, Greens ):
        """
        Procudes a matrix that maps from states where the excitation is matter
        to those where it is phtonic.
        (Opposite of fanout)
        """
        
        M_out, in_out, out_out = self.fanout(Greens)
        
        # Hermitian Conjugate of fanout, and swap in and out.
        return M_out.H, out_out, in_out
    
    
    def detuning_mat( self, inphoton_energy ):
        """
        Creates a diagonal square matrix of size equal to the number of
        excited states. Each element is the detuning between the given
        excited state and the initial state.
        
        If the system has ground states at multiple energies, then the 
        initial state energy is not fully determined by the input photon
        energy. Thus, a list of matrices are produced, one for each
        possibility.

        Parameters
        ----------
        inphoton_energy : float
            Energy of the input photon.

        Returns
        -------
        matrix_list : list of matrix
            List of detuning matrices.
        """
        
        matrix_list = []
        
        # list(set(...)) means I consider each energy once, even if it occurs
        # multiple times.
        for g_en in list(set(self.grd_ens)):
        
            det_list = []
            
            for e_en in self.ext_ens:
                det_list.append( e_en - inphoton_energy - g_en )
             
            matrix_list.append( np.matrix( np.diag( det_list ) ) )
                
        return matrix_list
    
    def ground_manifolds( self ):
        """
        Two ground states of different energies cannot be coupled by the long
        photon scattering. So need a filtering matrix to impose this.

        Returns
        -------
        filtering_matrix

        """
        n_gs = len(self.grd_ens)
        
        out_mat = np.zeros( (n_gs, n_gs) )
        
        for n1, en_g1 in enumerate(self.grd_ens):
            for n2, en_g2 in enumerate(self.grd_ens):
                if en_g1 == en_g2:
                    out_mat[n1, n2] = 1
                    
        return out_mat
    
    def ground_filters( self ):
        """
        Two ground states of different energies cannot be coupled by the long
        photon scattering. So need a filtering matrix to impose this.

        Returns
        -------
        filtering_matrix

        """
        n_gs = len(self.grd_ens)
        outlist = []
        for unique_energy in list(set(self.grd_ens)):
            outlist.append( [ item==unique_energy for item in self.grd_ens ] )
                    
        return outlist
    
    def inner(self, Greens_funct):
        """
        Calculates the coupling matrix between the excited states.
        Each element of this matrix gives the total (complex) coupling
        coeficient, this includes all possible ways for the first excited
        state to decay to a ground state with a photon emission, and for that
        ground state to re-absorb that photon and rise to the second excited
        state.

        Parameters
        ----------
        Greens_funct : Greens function object
            The Green's function describing the optical modes.

        Returns
        -------
        out_mat : matrix
            The combined coupling matrix.

        """
        out_mat = np.zeros((self.N_ext, self.N_ext), dtype= complex)
        
        for n in range(self.N_grd):
            # Gamma is the sum over W, taking the diagonals in the ground states.
            
            for G_n in range(len(Greens_funct.G_list)):
                # sum over THE DIAGONALS, IE the two ns are the same, and the two G_ns the same.
                out_mat += self.partial_inner( Greens_funct, n, n, G_n, G_n )
        
        return out_mat
    
    
    def partial_inner( self, Greens_funct, n1, n2, G_n1, G_n2 ):
        """
        Calculates a partial inner product.
        
        Produces a matrix, of size equal to the number of excited states square.
        Each matrix element represents the complex amplitude with which the
        first excited state can decay into ground state "n1", and mode "G_n1",
        and then for ground state "n2" with a photon in mode "G_n2" to excite
        the second excited state.
        Something like:
        <e2| H |g_n2, 1_{Gn2}> <g_n1, 1_{Gn1}| H |e1>
        
        Note, that it _is_ very weird that two different ground states 
        and two different optical modes are given, rather than repeating the 
        same one to make a valid path. Usually, this function _will_ be used
        with n1=n2 and G_n1=G_n2 giving this right-seeming behaviour.
        
        However, when calculating the "sandwhich" decay operators H\rhoH
        these non-diagonal terms do emerge.

        Parameters
        ----------
        Greens_funct : Greens function object
            The Green's function describing the optical modes.
        n1 : int
            The first ground state.
        n2 : int
            The second ground state.
        G_n1 : int
            The first optical mode.
        G_n2 : int
            the second optical mode.

        Returns
        -------
        matrix
            The coupling matrix.

        """
        
        # Like inner, but contrained to a particular pairing of ground states and modes.
        out_mat = np.zeros((self.N_ext, self.N_ext), dtype= complex)
        
        for i in range(self.N_ext):
            for j in range(self.N_ext):
                # For every pair of excited states i, j
                out_mat[i, j] = ( self.dipoles[i][n1] * 0.5 * Greens_funct.E_list[G_n1].H *  
                                    Greens_funct.E_list[G_n2] * self.dipoles[j][n2].H ).item()
        
        return np.matrix(out_mat)
    
    
    def represent(self):
        """Draws a cartoon of itself
        """
        
        fig, ax = plt.subplots(figsize = (5,5))
        xscale = max(self.ext_ens) - min(self.grd_ens)
        
        for x, (name, energy) in enumerate(zip(self.g_names, self.grd_ens)):
            ax.plot( [xscale*(x-0.2), xscale*(x+0.2)], [energy, energy], 'k-' )
            ax.text( xscale * x, energy - 0.1, name )
            
        for x, (name, energy) in enumerate(zip(self.e_names, self.ext_ens)):
            ax.plot( [xscale*(x-0.2), xscale*(x+0.2)], [energy, energy], 'k-' )
            ax.text( xscale * x, energy + 0.1, name )
            
        for i in range(len(self.e_names)):
            for j in range(len(self.g_names)):
                
                if np.linalg.norm( np.array( self.dipoles[i][j] )) > 0:
                    # It is nonzero!
                    #print(np.array(self.dipoles[i][j]))
                    #print(np.array(self.dipoles[i][j])[0, 0])
                    ax.plot( [xscale*i, xscale*j], [self.ext_ens[i], self.grd_ens[j]], 'r-' )
                    
                    draw_ellipse( np.array(self.dipoles[i][j])[0, 0], np.array(self.dipoles[i][j])[0, 1], 
                                 centre = [xscale * (i + j) / 2 , (self.ext_ens[i] + self.grd_ens[j]) / 2],
                                 scale = 0.8 * xscale, ax = ax )
        
        ax.set_aspect('equal', 'box')
        ax.axis('off')
        return ax


# def combine_systems( S1, S2 ):
#     # Given two atoms S1 and S2 I can construct a combined atom for the two of them.
#     # This is only valid in the single-excitition subspace.
#     import itertools
    
#     g_names = [ item1 +", " + item2 for item1, item2 in list(itertools.product( S1.g_names, S2.g_names) )]
#     e_names = [ item1 +", " + item2 for item1, item2 in list(itertools.product( S1.e_names, S2.g_names) )] + [ item1 +", " + item2 for item1, item2 in list(itertools.product( S1.g_names, S2.e_names) )]
    
#     # ground is them both in a ground state
#     g_ens = [ item1 + item2 for item1, item2 in list(itertools.product( S1.grd_ens, S2.grd_ens) )]
    
#     # excited is the ways of having just one excited.
#     e_ens = [ item1 + item2 for item1, item2 in list(itertools.product( S1.ext_ens, S2.grd_ens) )] + [ item1 + item2 for item1, item2 in list(itertools.product( S1.grd_ens, S2.ext_ens) )]
    
#     # Doubly excited manifold ignored by assumption, because we assume a single excitiaiton.
    
#     S_Out = multilevel_system( e_ens,
#                                g_ens,
#                                e_names = e_names, g_names = g_names ) # ground is the kron of both grounds.
    
#     ## Add transitions
#     for enum, e1 in enumerate(S1.e_names):
#         for gnum, g1 in enumerate(S1.g_names):
#             for g2 in S2.g_names:
#                 S_Out.add_transition( e_names.index(e1+", "+g2), g_names.index(g1+", "+g2), S1.dipoles[enum][gnum] )
    
#     for enum, e2 in enumerate(S2.e_names):
#         for gnum, g2 in enumerate(S2.g_names):
#             for g1 in S1.g_names:
#                 S_Out.add_transition( e_names.index(g1+", "+e2), g_names.index(g1+", "+g2), S2.dipoles[enum][gnum] )
                
#     return S_Out
    

# Note, that for scattering I am currently asusmign all ground states degenerate.
# I can do better than that.


def long_photon_scattering_matrix( G, S, freq = 1.0 ):
    """
    The scattering matrix of the atomic system S with the Greens function G.

    Parameters
    ----------
    G : Greens function object
        The Green's function describing the optical modes.
    S : Multilevel system object
        The Multilevel quantum system.
    freq : float
        The energy of the input photon.
    
    Returns
    -------
    scattering matrix : matrix
    inspace : str list,
        the state space
    """
    
    # Only start and end modes othe same energy are coupled in long-photon limit
    E_modes = len( G.G_list )
    
    Ground_energy_krondelta = np.kron( np.ones((E_modes,E_modes)), S.ground_manifolds() )
    
    Gamma = S.inner( G )
    
    fout, inspace, outspace = S.fanout( G )
    fin  = S.fanin( G )[0]
    
    # The detuning for each unique ground state energy.
    Det_list = S.detuning_mat( freq )
    filter_list = S.ground_filters()
    
    mat = np.eye( E_modes*len(S.grd_ens) , dtype = complex)
    
    
    for ground_filter, Det in zip(filter_list, Det_list):
        
        filtmat = np.matrix( np.diag( np.kron( np.ones((E_modes)), ground_filter) ) )
        
        # filtmat - fin * lng.inv( Gamma + 1j * Det ) * fout
        mat +=  -1 * fin * lng.pinv( Gamma + 1j * Det ) * fout * filtmat
        
    
    return mat, inspace


def long_photon_scattering( G, S, init_g, in_put, freq = 1.0 ):
    """
    Calculates the output state, after launching a photon at a system.

    Parameters
    ----------
    G : Greens function object
        The Green's function describing the optical modes.
    S : Multilevel system object
        The Multilevel quantum system.
    init_g : list
        List of initial complex amplitudes for the different ground states
    in_put : list
        List of initial complex amplitudes for the optical states
    freq : float
        The frequency of the input photon.
    
    Returns
    -------
    outstate : list of complex amplitudes
    state_space of output : list of strings
    """
    
    instate = np.matrix( np.kron( in_put, init_g ) )
    
    mat, inspace = long_photon_scattering_matrix( G, S, freq )
    
    return mat * instate.T, inspace
    
    
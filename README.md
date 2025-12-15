# GreenScat
Green'.s function based scattering model of waveguided photons from multilevel quantum systems.

The user provides information about a multilevel quantum system (for example an atom of a specific species), and the paramters of the optical modes the system is interacting with. The program then outputs scattering coefficients, assuming the system is prepared in a ground state and a single photon is thrown in (base use). The program is specifically built to deal with "wheel like" circular polarisation in the optical modes, which gives rise to interesting effects often called "chirality". 

If you also have the QuTiP package installed, the program can instead calculate the time evolution begining from an excited state.

## Example Use

You want to look at the situation described in doi:10.1126/science.1254699. Here, a Caesium atom interacts with a waveguide that posseses a circular polarisation (wheel-like circular polarisation, not bullet like circular as occurs in free space).

Import.

```
import MGS
```

Settup a multilevel system. As a simple representation of a Caesium atom we will assume a "Lambda"-like energy level configuration with one excited state and two ground states. We will set the energy of the excited state to be "10" units, and that of both ground states to "0".

```
S = MGS.multilevel_system([10], [0, 0])
```

Next, we tell the program which transitions between energy levels are dipole-allowed, and what dipoles they are associated with. In this case from excited state 0 to ground state 0 is allowed with a circular dipole, and from the excited state to the other ground state is allowed with the other helicty of circular dipole.

```
S.add_transition( 0, 0, [1, +1j] )
S.add_transition( 0, 1, [1, -1j] )
```

Next, we settup the optical enviroment of the atom. In this case it is a waveguide with a circular polarisation at the atom's location. (Note again, this is wheel-like circular polarisation). The interaction strength of the waveguide modes is determined largely by the group velocity (vg) that is assumed. We also include two loss modes, allowing the atom to loose light out of the waveguide. Two loss modes are included as this allows the loss to couple to all dipoles equally, which represents a simpler and more symmetric situation than would be possible with only one loss mode. The interaction strengths of the loss modes are set to 0.2 arbitary units.

```
vg = 0.1
mag = (0.5 / vg)
G = MGS.WG_Greens_function( [[1, 1j]], [mag], ["WG"])
G.add_mode( [1, 0], 0.2, False, "H loss" ) # Here "False" means it is not a waveguide mode. (Waveguide modes always come in pairs, one forward and one backward).
G.add_mode( [0, 1], 0.2, False, "V loss" )
```

Now we bounce a photon off of it. We will go for a photon resonant with the transitions (so 10 frequency/energy units). We will have the atom begin the interaction prpeared in its first ground state (indicated by the [1, 0], an arbitary initial superposition would be indicated [alpha, beta]). The photon enters from the forward mode ([1, 0, 0, 0], again an arbitary superpsotion of the forward, backward, and two loss modes is possible by replacing this with some vector of complex numbers). We will then print the resulting end state. 

```
scat, space = MGS.long_photon_scattering(G, S, [1, 0], [1, 0, 0, 0], freq=10)

print('Lambda System, circular dipoles, c point. Should toggle-reflect.')
for r, name in zip([complex(i[0,0]) for i in scat], space):
    print( name, " " * (24-len(name)), r)
```
This prints some numbers:

Lambda System, circular dipoles, c point. Should toggle-reflect.
WG-for, g1                (0.038461538461538325+0j)
WG-for, g2                0j
WG-back, g1               0j
WG-back, g2               (-0.9615384615384617+0j)
H loss, g1                (-0.1359820733051053+0j)
H loss, g2                (-0.1359820733051053+0j)
V loss, g1                -0.1359820733051053j
V loss, g2                0.1359820733051053j

We can see that the WG-back, g2 probability-amplitude reaches ~-0.96. Meaning a high probability that the photon is refelcted into the backward mode, with the atom moved simultaneously into the second ground state. This completes the simple example.



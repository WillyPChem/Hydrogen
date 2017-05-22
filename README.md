# Hydrogen
Simulator for Hydrogen in intense laser field.
Developed by Henrietta Eugene, Brynmor Anderson, Rana Daher, Laura Perdomo-Quintero, Muttalib Khan, Oyewumi Akinfenwa, and Jonathan J. Foley.

Uses time-dependent configuration interaction (TDCI) theory to solve the time-dependent Schrodinger equation for Hydrogen with dipolar coupling to an intense laser field.  The shape and intensity of the laser field can be adjusted, and the response of the hydrogen atom via its time-dependent dipole moment and the related absorption spectra can be computed.  The TDCI is performed in the basis of hydrogens energy eigenstates.  Recusrsive formulae are used for the computation of Laguerre polynomials in the radial component of the wavefunction, and for the Legendre polynomials in the angular part of the wavefunction.  A simple and relatively inefficient implementation of the discrete Fourier transform is used to compute the absorption spectrum from the time-dependent dipole moment. 

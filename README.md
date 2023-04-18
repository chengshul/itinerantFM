# iter_FM
Exact diagonalization (ED) and density-matrix renormalization group (DMRG) calculations of itinerant ferromagnetism in three systems.

## square_tri
The Hubbard model on the square lattice with unidirectional next-nearest-neighbor hopping, interpolating between the square lattice and the triagular lattice.

## square_flux
The Hubbard model on the square lattice with uniform commensurate flux.

## honeycomb
The Hubbard model on the honeycomb lattice with unidirectional next-nearest-neighbor hopping.

In all cases, ed\*.py calculates the case with two particles using ED, and dmrg\*.jl calculates the general filling with DMRG using the iTensor Library. In ED, we calculate the ground state energy in the singlet and triplet sectors. In DMRG, we calculate the ground state energy and total spin in the S_{z,tot}=0, S_{z,tot,max} sectors.

# itinerantFM
Exact diagonalization (ED) and density-matrix renormalization group (DMRG) calculations of itinerant ferromagnetism in three systems, arXiv:2305.01682.

## square_tri
The Hubbard model on the square lattice with unidirectional next-nearest-neighbor hopping, interpolating between the square lattice and the triagular lattice.

## square_flux
The Hubbard model on the square lattice with uniform commensurate fluxes.

## honeycomb
The Hubbard model on the honeycomb lattice with unidirectional next-nearest-neighbor hopping.

## codes 
In each case, the scripts `ed*.py` calculate the case with two particles using ED, and the scripts `dmrg*.jl` calculate the general filling with DMRG using the iTensor Library.
In ED, we calculate the ground state energy in the singlet and triplet sectors. In DMRG, we calculate the ground state energy and total spin in the $S_{z,tot}=0, S_{z,tot,max}$ sectors.

## data
In each case, DMRG data for energy and total spin are saved in `data/` folder.
The scripts `plot*.py` parse and plot the data of the given system size, particle number, etc. The scripts `plot_all.py` summarize and plot all data.

Packages
========

``LatticeQM`` provides several subpackages:

 * ``Structure`` contains several types and methods to define and
   manipulate lattice vectors, unit cells and paths in reciprocal space.
   It is an abstract class for many problems.
 * ``Geometries2D`` contains predefined two-dimensional lattice
   geometries (so far mostly honeycomb multilayers). It is a collection of
   tested examples.
 * ``TightBinding`` contains types and methods to define real-space hopping
   Hamiltonians (non-interacting), in particular via given ``Lattice`` objects.
   A tight-binding Hamiltonian can then be converted in to a Bloch Hamitlonian
   matrix. It is an abstract class intended to fit many problems.
 * ``Materials`` is a collection of physical examples and general modifiers to
   get tight-binding operators, e.g. for graphene, or zeeman fields, spin-orbit
   coupling etc.
 * ``BlochTools`` contains methods to perform, save and display results from
   exact diagonalization, such as bandstructure, optical conductivity,
   topological invariants, etc.

There is another subpackage ``KPM.jl`` implementing the Kernel Polynomial Method.
It is tested and working, but for now neither integrated nor properly documented.

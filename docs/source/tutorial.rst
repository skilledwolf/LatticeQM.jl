Tutorial
========

The typical workflow is as follows:

1. Import the package via ``using LatticeQM``

2. Create or load a predefined lattice object

3. Obtain a tight-binding Hamiltonian from the lattice object

4. Use exact diagonalization to calculate bandstructure,
   linear response coefficients and topological invariants

In what follows, we shall elaborate each point. Afterwards, you will have
a good idea of how to use this package.

Defining a lattice
------------------

A :math:`d`-dimensional lattice in :math:`(d+D)`-dimensional space is specified through its lattice vectors
:math:`\vec{a}_i\in\mathbb{R}^{d+D}` for :math:`i=1,\dots,d`.
We further need to specify the atom (and/or orbital) positions :math:`\vec{r}_\alpha=\sum_i (\vec{x}_\alpha)_i \vec{a}_i` for :math:`\alpha=1,\dots,N`
per unit cell, in terms of fractional coordinates :math:`\vec{x}_\alpha`.

Lattice examples
""""""""""""""""
As example, we define a two-dimensional two-atomic square lattice with explicit
z-coordinates via

.. code-block:: julia
  :linenos:

  extra_dimensions = ["z"] # allows to specify additional coordinates

  A = [[1, 0, 0] [0, 1, 0]] # each column is a lattice vector
  atoms = [[0, 0, 0] [0.5, 0.5, 1.0]] # each column is an atom within the unit cell (in fractional coordinates)

  lat = Lattice(A, atoms, extra_dimensions)

The case without extra coordinates is

.. code-block:: julia
  :linenos:

  A = [[1, 0] [0, 1]] # each column is a lattice vector
  atoms = [[0, 0] [0.5, 0.5]] # each column is an atom within the unit cell (in fractional coordinates)

  lat = Lattice(A, atoms)

The case with a single atom per unit cell is

.. code-block:: julia
  :linenos:

  A = [[1, 0] [0, 1]] # each column is a lattice vector
  lat = Lattice(A)

or even shorter

.. code-block:: julia
  :linenos:

  lat = Lattice(2)

Modifying the ``Lattice`` object after creation is discouraged.
Use the methods ``update_atoms!(lat, atoms)`` and ``update_A!(lat,A)`` if you
really have to.

To access lattice properties use ``get_A(lat)``, ``get_coordinates(lat)``,
``get_positions(lat)``, ``get_positions_in(lat, extradim_name)``,
``get_extrapositions(lat)``, ``get_ldim(lat)`` and ``get_sdim(lat)``.
Reciprocal lattice vectors can be obtained as columns of ``get_B(lat)``.

Several predefined lattices can be found in ``LatticeQM.Geometries2D``, e.g.,

.. code-block:: julia

  lat1 = Geometries2D.honeycomb()
  lat2 = Geometries2D.honeycomb_twisted(7)
  lat3 = Geometries2D.honeycomb_twisted_ABAB(7)

Visualizing a two-dimensional lattice
"""""""""""""""""""""""""""""""""""""
Currently only two-dimensional lattices can be visualized. In the future, it
would be nice if higher-dimensional lattices could be visualized via 2D
projections.

.. code-block:: julia

  lat = Geometries2D.honeycomb_twisted(7)
  plot(lat; repeat=[0:1,0:1])


Creating a tight-binding Hamiltonian
------------------------------------

A translationally-invariant Hamiltonian is given by

.. math::

	H = \sum_{\vec{r},\delta\vec{r}\in\Lambda_d} \sum_{\alpha,\beta=1,\dots,N} \sum_{s,s'} t_{\delta\vec{r}}^{\alpha\beta, s s'} \; c_{\vec{r}+\delta\vec{r},\alpha s}^\dagger c_{\vec{r},\beta s'} ,

where :math:`\Lambda_d` is the set of lattice points and
:math:`c_{\vec{r},\alpha s}^{(\dagger)}` destroys (creates) an electron with
spin :math:`s` at orbital :math:`\alpha` in unit cell :math:`\vec{r}`.

Clearly, the Hamiltonian is fully specified by the hopping matrices
:math:`t_{\delta\vec{r}}`. Typically, only short-range hoppings are relevant.

Via Lattice object and distance function
""""""""""""""""""""""""""""""""""""""""

The most convenient way to define a tight-binding Hamiltonian is via the
``Lattice`` object and a distance function :math:`f(\vec{r}_1,\vec{r}_2)` that
returns the hopping amplitude between atoms/orbitals at positions
:math:`\vec{r}_1` and :math:`\vec{r}_2`.

In the following example, we set nearest-neighbor hopping amplitudes to 1, while
all others are zero:

.. code-block:: julia
  :linenos:

  lat = Lattice(2) # single atom square lattice

  f(r1,r2=0.0) = (dot(r1.-r2,r1.-r2)<1.1) ? 1.0 : 0.0 # hopping function

  hops = get_hops(lat, f) # get hopping matrix

Explicit definition
"""""""""""""""""""

For small problems, it might be ok to set the hopping elements by hand.

For example, the tight-binding Hamiltonian for a minimal model of graphene is

.. code-block:: julia
  :linenos:

  t0 = complex([ 0.0  1.0;
                 0.0  0.0 ])

  hops = Dict(
    [0, 0] => t0 + t0',
    [0, 1] => t0,
    [1, 0] => t0,
    [0, -1] => t0'
    [-1, 0] => t0'
  )

Predefined tight-binding terms
""""""""""""""""""""""""""""""

.. code-block:: julia

  <code>


Bandstructure calculations
--------------------------

Text here.

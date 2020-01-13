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

To access lattice properties use ``getA(lat)``, ``coordinates(lat)``,
``get_positions(lat)``, ``extrapositions(lat, extradim_name)``,
``get_extrapositions(lat)``, ``get_ldim(lat)`` and ``get_sdim(lat)``.
Reciprocal lattice vectors can be obtained as columns of ``getB(lat)``.

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

  f(r1,r2=0.0) = (dot(r1.-r2,r1.-r2)<1.1) ? 1.0 : 0.0 # hopping function

  lat = Lattice(2) # single atom square lattice
  hops = gethops(lat, f) # get hopping matrix

or

.. code-block:: julia
  :linenos:

  f(r1,r2=0.0) = (dot(r1.-r2,r1.-r2)<1.1) ? 1.0 : 0.0 # hopping function

  lat = Lattice(2)
  hops = Hops()

  addhops!(hops, lat, f)



Explicit definition
"""""""""""""""""""

For small problems, it might be ok to set the hopping elements by hand.

For example, the tight-binding Hamiltonian for a minimal model of graphene is

.. code-block:: julia
  :linenos:

  t0 = [ 0.0  1.0;
         0.0  0.0 ]

  hops = Hops(
    [0, 0] => t0 + t0',
    [0, 1] => t0,
    [1, 0] => t0,
    [0, -1] => t0'
  )

  hops[[-1,0]] = t0' # hops can be set/modified any time


Predefined tight-binding terms
""""""""""""""""""""""""""""""

There are several predefined tight-binding terms, such as:

.. code-block:: julia

  hops = Materials.graphene(lat; mode=:spinhalf)

  TightBinding.add_chemicalpotential!(hops, lat, 0.1)
  TightBinding.set_filling!(hops, lat, 0.5; nk=100)

  Materials.addsublatticeimbalance!(hops, lat, 0.1)
  Materials.addhaldane!(hops, lat, 0.2; spinhalf=true)
  Materials.add_spinorbit!(hops, lat, 0.03)
  Materials.addrashba!(hops, lat, 0.2)
  Materials.add_transversepotential!(hops, lat, 0.1)
  Materials.add_zeeman!(hops, lat, [0.0,0.0,0.1])


Calculations
------------

Bandstructure
"""""""""""""

.. code-block:: julia

  lat = Geometries2D.honeycomb()
  hops = Materials.graphene(lat; mode=:spinhalf)
  h = getbloch(hops)

  # Build up a valley-observable
  valleyhops = Materials.get_haldane_hops(lat, √3/9; spinhalf=true, mode=:sublatticeA)
  Materials.addhaldane!(valleyhops, lat, -√3/9; spinhalf=true, mode=:sublatticeB)
  valley_proj = TightBinding.expvalf(getbloch(valleyhops))

  # ks = Structure.kpath(lat; num_points=200)
  ks = kpath(Geometries2D.k_hexagonal; num_points=200)
  bands = get_bands(h, ks; projector=valley_proj) #getprojector(lat, "spin")

  # Show bands
  save(bands, "graphene_bands.h5")
  plot(bands, ylabel="\$\\varepsilon/t\$", colorbar_title="valley", size=(330,240), colorbar=true, markercolor=:PiYG)


Density of states
"""""""""""""""""

.. code-block:: julia

  nk = 400^2
  N = 500
  Λ = 3.1

  ks = Structure.regulargrid(nk=nk)
  energies = collect(range(-Λ, length=N, stop=Λ))
  Γ = 0.01
  DOS = BlochTools.dos_dense(h, ks, energies; Γ=Γ)

  plot(energies, DOS)

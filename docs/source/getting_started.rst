Getting started
===============

Installation
------------

 * Make sure `Julia <https://julialang.org>`_ is installed and is at least version 1.2.

 * Use Julia's builtin package manager `Pkg` to add
   `LatticeQM.jl <https://gitlab.ethz.ch/wolft/LatticeQM.jl>`_. Note that at the
   time of writing, the repository is private and requires login credentials.
   Just contact me if you want access.

   **Option 1:** Start an interactive julia session (`REPL`), hit the ``]`` key
   and run::

      add https://gitlab.ethz.ch/wolft/LatticeQM.jl

   Once the installation is done, leave Pkg mode with the `backspace` key.

   **Option 2:** Execute the julia code

   .. code-block:: julia

      using Pkg
      Pkg.add("https://gitlab.ethz.ch/wolft/LatticeQM.jl.git")

   **Option 3:** If you want to develop the package, you could do 

   .. code-block:: bash

      $ cd /where/I/want/my/files/
      $ git clone https://gitlab.ethz.ch/wolft/LatticeQM.jl.git LatticeQM
      $ julia
      julia> ]
      pkg> develop /where/I/want/my/files/LatticeQM


Usage
-----

.. code-block:: julia

  using LatticeQM

Note that the first import may take quite long (due to compilation). If you get
error messages about third-party packages try to install them manually through
the Pkg package manager. If it ran successfully, then you are good to go.

Some functions, modules and types are documented with docstrings, in which case you can
run

.. code-block:: bash

   julia> ?getbands

which will display a description of the function and its arguments together with a short
usage example.

First-time users  may want to check out the tutorials and examples.

Example code
------------

.. literalinclude:: codes/example.jl
   :language: julia

The output is

.. image:: codes/example.png
  :width: 400
  :alt: Band structure for graphene
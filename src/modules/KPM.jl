# __precompile__()

module KPM

    """
    This module implements the "Kernel Polynomial Method" to obtain expectation
    values and stochastic traces traces of Hermitian matrices with bounded spectra.

    This method is particularly suitable for sparse 0D systems. It does not
    require diagonalization (except to find largest and smallest eigenvalue).
    The computational cost increases only linearly in system size and the
    expansion is numerically stable. The stochastic trace requires fewer random
    vectors for larger systems.

    For details on how the method works and for which problems it can be applied,
    please refer to the review paper Ref. [1].

    [1] Weisse et al, Rev. Mod. Phys. 78 275
    """

    include("KPM/misc.jl")
    include("KPM/chebyshev.jl")
    include("KPM/kernels.jl")
    include("KPM/recursion.jl")

    include("KPM/applications.jl")

end

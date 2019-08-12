using Printf

function search_fixedpoint!(f!, x1, x0;
    iterations=500,
    tol=1e-7, β=1.0,
    p_norm::Real=2,
    show_trace=false)
    # show_report=false)
    """
        Fixedpoint iteration, tested on the square-root example
        f_a(x) = 1/2 * (a/x+x)
        which has the fixed point x0 = sqrt(a).
    """

    converged = false
    ϵ0 = 0.0

    if show_trace #|| show_report
        println("==============================")
        println(" FIXPOINT SEARCH ")
        println(" #  \t error")
        println("==============================")
    end

    ϵ_abs = 1.0
    iter = 0
    while iter < iterations
        iter += 1

        # Perform a step
        ϵ0 = f!(x1, x0)

        # Convergence?
        ϵ_abs = norm(values(x1).-values(x0), p_norm)

        if show_trace
            print("\r")
            print(@sprintf(" %d  \t %.2E", iter, ϵ_abs))
            print("\u1b[0K")
        end

        if ϵ_abs < tol
            converged = true
            break
        end

        # Convergence acceleration for the next step
        # The new x0 for the next step is:
        for δL=keys(x0)
            @. x0[δL] .= β * x1[δL] + (1-β) * x0[δL]
        end
        # @. x0 .= (1-β) * x0 + β * x1
    end


    # if show_report && !show_trace
    #     println(@sprintf(" %d  \t %.2E", iter, ϵ_abs))
    # end

    if show_trace #|| show_report
        if converged
            println("\nConverged!\n")
        else
            println("\nNOT converged.\n")
        end
    end

    ϵ0, converged
end

################################################################################
################################################################################
################################################################################
################################################################################

function solve_selfconsistent(ℋ_op::Function, ℋ_scalar::Function,
    ρ_init::Dict{Vector{Int},AbstractArray{T0,N}}, ks::AbstractMatrix{Float64}, filling::Float64;
    iterations=500, tol=1e-7, T=0.0, format=:dense, kwargs...) where {N, T0<:Number}
    """
        Searches a self-consistent meanfield solution for the functional

            ℋ: ρ → h  where h(k) is a hermitian N × N Matrix

        at given filling (between 0 and 1). k space is discretized with
        the given points ks (columns).

        Note: this amounts to a fixed-point search.
    """

    type = issparse(ℋ_op(ρ_init)(ks[:,1])) ? :sparse : :dense # Decide if the Hamiltonian is sparse

    function update_ρ!(ρ1, ρ0)

        # Update meanfield Hamiltonian and chemical potential
        h = ℋ_op(ρ0) # probably o.k.

        Σ = spectrum(h; format=format)
        μ = chemical_potential(h, ks, filling; T=T)#; type=type)

        # Obtain the meanfield density matrix of the updated Hamiltonian
        ϵ0 = ρ_L!(ρ1, Σ, ks, μ; T=T)

        ϵ0 # return the groundstate energy (density matrix was written to ρ1)
    end

    # Compute the ground state energy for the mean-field fixed point
    ρ0 = deepcopy(ρ_init)
    ρ1 = deepcopy(ρ_init)

    ϵ0, converged = search_fixedpoint!(update_ρ!, ρ1, ρ0; iterations=iterations, tol=tol, kwargs...)

    ρBloch = build_BlochH(ρ1; mode=:nospin)

    hBloch = ℋ_op(ρ1)
    ϵ_offset = ℋ_scalar(ρ1)

    ϵ_GS = ϵ0 + ϵ_offset

    ρ1, ρBloch, hBloch, ϵ_GS, converged
end


################################################################################
################################################################################
################################################################################
################################################################################

function UnitMatrix(n::Int, i::Int)
    mat = spzeros(ComplexF64, n, n)
    mat[i,i] = 1.0+0.0im

    mat
end

function initialize_ρ(v::Dict{Vector{Int},T2}, mode=:random; lat=:nothing) where {T1<:Complex, T2<:AbstractMatrix{T1}}

    N = size(first(values(v)), 1)

    ρs = Dict{Vector{Int},AbstractMatrix{ComplexF64}}()
    for δL=keys(v)
        ρs[δL] = zeros(ComplexF64, size(v[δL]))
    end

    if mode==:randombig
        mat = rand(ComplexF64, N, N)
        ρs[zero(first(keys(ρs)))] = (mat + mat') ./ 2

    elseif mode==:random
        @assert mod(N,2)==0
        n = div(N,2)
        mat = rand(ComplexF64, n, n)

        # Generate a random spin orientation at a lattice site
        function randmat()
            # d = rand(Float64, 3)
            d = -1.0 .+ 2 .* rand(Float64, 3)
            p = 0.5 .* (σ0 .+ sum(d[i_]/norm(d) .* σs[i_] for i_=1:3))
        end

        ρs[zero(first(keys(ρs)))] = Matrix(sum(kron(randmat(),UnitMatrix(n,i_)) for i_=1:n))

    elseif mode==:antiferro || mode==:antiferroZ
        sublA, sublB = get_operator(lat, ["sublatticeA", "sublatticeB"])
        mat = sublA .- sublB

        σUP = 0.5 .* (σ0 .+ σZ)

        ρs[zero(first(keys(ρs)))] = kron(σUP, mat)

    elseif mode==:antiferroX
        sublA, sublB = get_operator(lat, ["sublatticeA", "sublatticeB"])
        mat = sublA .- sublB

        σUP = 0.5 .* (σ0 .+ σX)

        ρs[zero(first(keys(ρs)))] = kron(σUP, mat)

    elseif mode==:ferro || mode==:ferroZ #|| mode==:ferroz
        @assert mod(N,2)==0
        n = div(N,2)
        σUP = 0.5 .* (σ0 .+ σZ)
        ρs[zero(first(keys(ρs)))] =  2. * kron(σUP, Diagonal(ones(n)))

    elseif mode==:ferroX #|| mode==:ferroz
        @assert mod(N,2)==0
        n = div(N,2)
        σLEFT = 0.5 .* (σ0 .+ σX)
        ρs[zero(first(keys(ρs)))] =  2. * kron(σLEFT, Diagonal(ones(n)))
    # elseif mode==:ferrox
    #     @assert mod(N,2)==0
    #     n = div(N,2)
    #     ρs[zero(first(keys(ρs)))] = kron(σX, Diagonal(ones(n)))

    else
        error("Unrecognized mode '$mode' in initialize_ρ(...).")

    end

    ρs
end



################################################################################
################################################################################
################################################################################
################################################################################

getdiag(A::AbstractMatrix) = view(A,diagind(A,0))
BlochPhase(k,δL)::ComplexF64  = exp(1.0im * 2 * π * ComplexF64(dot(k,δL)))

function get_mf_functional(h::Function, v::Dict{Vector{Int},T2}; format=:sparse) where {T1<:Complex, T2<:AbstractMatrix{T1}}
    """
        This method takes the Hamiltonian single-particle operator h and an
        interaction potential v and returns mean-field functionals
            ℋ, E s.t. h_mf = ℋ[ρ] and ϵ_gs = E[ρ].

        These functionals can be used to search for a self-consistent solution
        using solve_selfconsistent(...).
    """

    mf_op, E = get_mf_operator(v; format=format)
    ℋ(ρ) = k -> (h(k) .+ mf_op(ρ, k))

    ℋ, E
end

function get_mf_operator(v::Dict{Vector{Int},T2}; format=:sparse) where {T1<:Complex, T2<:AbstractMatrix{T1}}
    """
        Expects the real space potential {V(L) | L unit cell vector}.
        It returns a functional 𝒱[ρ,k] that builds the mean field hamiltonian
        (i.e. h_v(k) = 𝒱[ρ,k]).

        This may look harmless but requires a careful derivation.
    """

    # d = size(first(values(v)),1)
    # vsym(L::Vector{Int}) = 0.5 .* (v[L].+(v[L])')
    V0 = sum(v[L] for L in keys(v))

    diag0(ρs) = getdiag(ρs[[0,0]])

    function mf_op(ρs::Dict{Vector{Int},T2}, k::AbstractVector{Float64}) where {T1<:Complex, T2<:AbstractMatrix{T1}}

        # Hartree contribution
        H_hartree = spdiagm(0 => V0 * (diag0(ρs))) #.-1/2

        # Fock contribution
        H_fock(k) = - sum(v[L] .* ρL .* BlochPhase(k,L) for (L,ρL) in ρs)

        H_hartree + H_fock(k) #+ real(e_hartree) * I + real(e_fock) * I
    end

    function mf_scalar(ρs::Dict{Vector{Int},T2}) where {T1<:Complex, T2<:AbstractMatrix{T1}}

        # Hartree contribution
        vρ = diag0(ρs)
        e_hartree = - 1/2 * (vρ' * V0 * vρ)
        @assert imag(e_hartree) ≈ 0

        # Fock contribution
        e_fock =  1/2 * sum(sum(ρL .* conj.(ρL) .* v[L] for (L,ρL) in ρs))
        @assert imag(e_hartree) ≈ 0

        real(e_hartree + e_fock)
    end

    mf_op, mf_scalar
end

# function get_mf_operator(v::Dict{Vector{Int},T2}; format=:sparse) where {T1<:Complex, T2<:AbstractMatrix{T1}}
#     """
#         Expects the real space potential {V(L) | L unit cell vector}.
#         It returns a functional 𝒱[ρ,k] that builds the mean field hamiltonian
#         (i.e. h_v(k) = 𝒱[ρ,k]).
#
#         This may look harmless but requires a careful derivation.
#     """
#
#     d = size(first(values(v)),1)
#     vsym(L::Vector{Int}) = 0.5 .* (v[L].+(v[L])')
#
#     function mf_op(ρs::Dict{Vector{Int},T2}, k::AbstractVector{Float64}) where {T1<:Complex, T2<:AbstractMatrix{T1}}
#
#         # Hartree contribution
#         vρ = getdiag(ρs[[0,0]])
#
#         H_hartree = spdiagm(0 => vsym([0,0]) * vρ)
#         # e_hartree = - 1/2 * (vρ' * v[[0,0]] * vρ)
#         # @assert imag(e_hartree) ≈ 0
#
#         # Fock contribution
#         H_fock(k) = - sum(vsym(L) .* transpose(ρL) .* BlochPhase(k,L) for (L,ρL) in ρs)
#         #e_fock =  1/2 * sum(sum(ρL .* transpose(ρL) .* vsym(L) for (L,ρL) in ρs))
#         # @assert imag(e_hartree) ≈ 0
#
#         H_hartree + H_fock(k) #+ real(e_hartree) * I + real(e_fock) * I
#     end
#
#     function mf_scalar(ρs::Dict{Vector{Int},T2}) where {T1<:Complex, T2<:AbstractMatrix{T1}}
#
#         # Hartree contribution
#         vρ = getdiag(ρs[[0,0]])
#
#         e_hartree = - 1/2 * (vρ' * v[[0,0]] * vρ)
#         @assert imag(e_hartree) ≈ 0
#
#         # Fock contribution
#         e_fock =  1/2 * sum(sum(ρL .* transpose(ρL) .* v[L] for (L,ρL) in ρs))
#         # e_fock = 0.0
#         @assert imag(e_hartree) ≈ 0
#
#         real(e_hartree + e_fock)
#     end
#
#     mf_op, mf_scalar
# end

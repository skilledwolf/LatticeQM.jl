using Printf

function search_fixedpoint!(f!, x1, x0;
    iterations=500,
    tol=1e-7, β=1.0,
    p_norm::Real=2,
    show_trace=false,
    clear_trace=false
    )
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
        println(" #  \t error \t time/step [s]")
        println("==============================")
    end

    error = 1.0
    iter = 0
    while iter < iterations
        iter += 1

        # Perform a step
        t0 = time_ns()
#         timediff = @timed ϵ0 = f!(x1, x0)
        ϵ0 = f!(x1, x0)
        t1 = (time_ns()-t0)/1e9

        # Convergence?
        error = norm(values(x1).-values(x0), p_norm)

        if show_trace
            if clear_trace print("\r") end
            print(@sprintf(" %d  \t %.2E \t %.2E", iter, error, t1))
            if clear_trace print("\u1b[0K") else println("") end
        end

        if error < tol
            converged = true
            break
        end

        if iter == iterations
            break
        end

        # Convergence acceleration for the next step
        # The new x0 for the next step is:
        for δL=keys(x0)
            @. x0[δL] .= β * x1[δL] + (1-β) * x0[δL]
        end
    end

    if show_trace #|| show_report
        if converged
            println("\nConverged!\n")
        else
            println("\nNOT converged.\n")
        end
    end

    ϵ0, error, converged
end

################################################################################
################################################################################
################################################################################
################################################################################
using SharedArrays

function solve_selfconsistent(ℋ_op::Function, ℋ_scalar::Function,
    ρ_init::Dict{Vector{Int},T1}, ks::AbstractMatrix{Float64}, filling::Float64;
    iterations=500, tol=1e-7, T=0.0, format=:dense, verbose::Bool=false, kwargs...) where {N, T0<:Number, T1<:AbstractArray{T0,N}}
    """
        Searches a self-consistent meanfield solution for the functional

            ℋ: ρ → h  where h(k) is a hermitian N × N Matrix

        at given filling (between 0 and 1). k space is discretized with
        the given points ks (columns).

        returns
            1) the density matrix of the meanfield
            2) ground state energy of the meanfield operator
            3) the chemical potential
            3) convergence flag (bool)
            4) error estimate

        Note: this amounts to a fixed-point search.
    """

    function update_ρ!(ρ1, ρ0)

        # Update meanfield Hamiltonian and chemical potential
        h = ℋ_op(ρ0)
        Σ = spectrum(h; format=:dense) # lazy diagonalization

        if verbose
            @info("Updating chemical potential for given filling.")
        end
        μ = chemical_potential(h, ks, filling; T=T)#; type=type) # @time

        # Obtain the meanfield density matrix of the updated Hamiltonian
        if verbose
            @info("Updating the meanfield density matrix.")
        end
        ϵ0 = ρ_L!(ρ1, Σ, ks, μ; T=T) # @time

        ϵ0 # return the groundstate energy (density matrix was written to ρ1)
    end

    # Compute the ground state energy for the mean-field fixed point
    ρ0 = Dict(δL=>SharedArray(m) for (δL, m)=ρ_init) #deepcopy(ρ_init)
    ρ1 = Dict(δL=>SharedArray(m) for (δL, m)=ρ_init)

    ϵ0, error, converged = search_fixedpoint!(update_ρ!, ρ1, ρ0; iterations=iterations, tol=tol, kwargs...)

    h = ℋ_op(ρ1)
    μ = chemical_potential(h, ks, filling; T=T) # Calculate the chemical potential at the end of iteration

    ϵ_offset = ℋ_scalar(ρ1)
    ϵ_GS = ϵ0 + ϵ_offset

    ρ1, ϵ_GS, μ, converged, error
end


################################################################################
################################################################################
################################################################################
################################################################################

BlochPhase(k,δL)::ComplexF64  = exp(1.0im * 2 * π * ComplexF64(dot(k,δL)))

function get_mf_functional(h::Function, v::Dict{Vector{Int},<:AbstractMatrix})
    """
        This method takes the Hamiltonian single-particle operator h and an
        interaction potential v and returns mean-field functionals
            ℋ, E  s.t.  h_mf = ℋ[ρ]  and  ϵ_scalar = E[ρ].

        These functionals can be used to search for a self-consistent solution
        using solve_selfconsistent(...).
    """

    mf_op, E = get_mf_operator(v)
    ℋ(ρ) = k -> (h(k) .+ mf_op(ρ, k))

    ℋ, E
end

function get_mf_operator(v::Dict{Vector{Int},<:AbstractMatrix})
    """
        Expects the real space potential {V(L) | L unit cell vector}.
        It returns a functional 𝒱[ρ,k] that builds the mean field hamiltonian
        (i.e. h_v(k) = 𝒱[ρ,k]).

        This may look harmless but requires a careful derivation.
    """

    V0 = sum(v[L] for L in keys(v))
    diag0(ρs) = diag(ρs[[0,0]])

    function mf_op(ρs::Dict{Vector{Int},<:AbstractMatrix}, k::AbstractVector{Float64})

        # Hartree contribution
        H_hartree = spdiagm(0 => V0 * (diag0(ρs)))

        # Fock contribution
        H_fock(k) = - sum(v[L] .* ρL .* BlochPhase(k,L) for (L,ρL) in ρs)

        H_hartree + H_fock(k)
    end

    function mf_scalar(ρs::Dict{Vector{Int},<:AbstractMatrix})

        # Hartree contribution
        vρ = diag0(ρs)
        e_hartree = - 1/2 * (transpose(vρ) * V0 * vρ)
        @assert imag(e_hartree) ≈ 0

        # Fock contribution
        e_fock =  1/2 * sum(sum(ρL .* conj.(ρL) .* v[L] for (L,ρL) in ρs))
        @assert imag(e_hartree) ≈ 0

        real(e_hartree + e_fock)
    end

    mf_op, mf_scalar
end

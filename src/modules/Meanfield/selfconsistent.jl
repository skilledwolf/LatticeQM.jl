
function solvemf_hubbard_k(lat, args...; kwargs...)

    v = gethubbard(lat)

    ℋ_op, ℋ_scalar = hartreefock(h, v)

    solveselfconsistent(ρ_init, ℋ_op, ℋ_scalar, ks, filling)
end

function solvemf_hartreefock_k(v::AnyHops, h::Function, ρ_init::AnyHops, ks::AbstractMatrix, filling::Number; kwargs...)
    ℋ_op, ℋ_scalar = hartreefock(h, v)

    solveselfconsistent(ρ_init, ℋ_op, ℋ_scalar, ks, filling; kwargs...)
end

solveselfconsistent(ρ_init, args...; kwargs...) = solveselfconsistent!(copy(ρ_init), args...; kwargs...)
solveselfconsistent!(ρ0::Dict{Vector{Int},T1}, ℋ_op::Function, ℋ_scalar::Function, ks::AbstractMatrix{Float64}, filling::Float64; kwargs...) where {T1<:AbstractArray} = solveselfconsistent!(ρ0, copy(ρ0), ℋ_op, ℋ_scalar, ks, filling; kwargs...)
function solveselfconsistent!(ρ0::Dict{Vector{Int},T1}, ρ1::Dict{Vector{Int},T1}, ℋ_op::Function, ℋ_scalar::Function,
    ks::AbstractMatrix{Float64}, filling::Float64;
    iterations=500, tol=1e-7, T=0.0, format=:dense, verbose::Bool=false, kwargs...) where {T1<:AbstractArray}
    """
        Searches a self-consistent meanfield solution for the functional

            ℋ: ρ → h  where h(k) is a hermitian N × N Matrix

        at given filling (between 0 and 1). k space is discretized with the given points ks.

        returns (1) the density matrix of the meanfield (2) ground state energy of the meanfield operator (3) the chemical potential (4) convergence flag (bool) (5) error estimate
    """
    ρ0 = Dict(δL=>SharedArray(m) for (δL, m)=ρ0)
    ρ1 = Dict(δL=>SharedArray(m) for (δL, m)=ρ1)

    function update!(ρ1, ρ0)
        # Update meanfield Hamiltonian and chemical potential
        h = ℋ_op(ρ0)
        Σ = spectrum(h; format=:dense) # lazy diagonalization

        if verbose
            @info("Updating chemical potential for given filling.")
        end
        μ = chemicalpotential(h, ks, filling; T=T)#; type=type) # @time

        # Obtain the meanfield density matrix of the updated Hamiltonian
        if verbose
            @info("Updating the meanfield density matrix.")
        end
        ϵ0 = densitymatrix!(ρ1, Σ, ks, μ; T=T) # @time

        ϵ0 # return the groundstate energy (density matrix was written to ρ1)
    end

    # Compute the ground state energy for the mean-field fixed point
    ϵ0, error, converged = fixedpoint!(update!, ρ1, ρ0; iterations=iterations, tol=tol, kwargs...)

    h = ℋ_op(ρ1)
    μ = chemicalpotential(h, ks, filling; T=T) # Calculate the chemical potential at the end of iteration

    ϵ_offset = ℋ_scalar(ρ1)
    ϵ_GS = ϵ0 + ϵ_offset

    ρ1, ϵ_GS, μ, converged, error
end

###################################################################################################
# Backwards compatibility
###################################################################################################
export solve_selfconsistent
@legacyalias solveselfconsistent solve_selfconsistent

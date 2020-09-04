using ..Utils: regulargrid

solvehartreefockmf_k(v::AnyHops, hops::AnyHops, args...; kwargs...) = solvehartreefockmf_k(v, getbloch(hops), args...; kwargs...)
function solvehartreefockmf_k(v::AnyHops, h::Function, ρ_init::AnyHops, filling::Number, args...; kwargs...)
    ℋ_op, ℋ_scalar = hartreefock(h, v)

    solveselfconsistent(ρ_init, ℋ_op, ℋ_scalar, filling, args...; kwargs...)
end

function solveselfconsistent(hf, ρ_init::AnyHops, filling::Number, args...; kwargs...)
    return solveselfconsistent(ρ_init, hf..., filling, args...; kwargs...)
end

function solveselfconsistent(ρ_init::AnyHops, args...; kwargs...)
    return solveselfconsistent!(deepcopy(ρ_init), args...; kwargs...)
end

function solveselfconsistent!(ρ0::AnyHops, ℋ_op::Function, ℋ_scalar::Function, filling::Float64, args...; kwargs...)
    return solveselfconsistent!(ρ0, deepcopy(ρ0), ℋ_op, ℋ_scalar, filling, args...; kwargs...)
end

function solveselfconsistent!(ρ0::AnyHops, ρ1::AnyHops, ℋ_op::Function, ℋ_scalar::Function, filling::Float64; klin::Int, kwargs...)
    solveselfconsistent!(ρ0, ρ1, ℋ_op, ℋ_scalar, filling, regulargrid(nk=klin^2); kwargs...)
end

mutable struct Hamiltonian
    h
    μ::Float64
end


using JLD


"""
    Searches a self-consistent meanfield solution for the functional ℋ: ρ → h
    at given filling (between 0 and 1). k space is discretized with the given points ks.
    returns (1) the density matrix of the meanfield (2) ground state energy of the meanfield operator (3) the chemical potential (4) convergence flag (bool) (5) error estimate

    If the checkpoint keyword is set, e.g. `checkpoint="mf.jld"`, the mean field will be saved into a file at each step of the iteration,
    allowing to interrupt and restart a long-running calculation safely. If `checkpoint` is set and the file exists, this method
    will assume that the file contains a valid mean-field and use it as initial guess. To avoid this behavior, specify `hotstart=false`.

    parallel=true might help if diagonalization per k point is very time consuming
    (e.g. for twisted bilayer graphene)
    note that for small problems `parallel=true` may decrease performance (communication overhead)
"""
function solveselfconsistent!(ρ0::AnyHops, ρ1::AnyHops, ℋ_op::Function, ℋ_scalar::Function, filling::Float64, ks::AbstractMatrix{Float64};
    convergenceerror=false, multimode=:serial, checkpoint::String="", hotstart=true, iterations=500, tol=1e-7, T=0.0, format=:dense, verbose::Bool=false, kwargs...)

    if checkpoint != "" && isfile(checkpoint) && hotstart
        println("Loading checkpoint file as initial guess: $checkpoint")
        ρ0 = JLD.load(checkpoint, "mf")
    end

    if multimode==:parallel
        ρ0 = Hops(δL=>SharedArray(complex(m)) for (δL, m)=ρ0) # ensure dense
        ρ1 = Hops(δL=>SharedArray(complex(m)) for (δL, m)=ρ1) # convert to dense
    else
        ρ0 = Hops(δL=>Matrix(complex(m)) for (δL, m)=ρ0) # ensure dense
        ρ1 = Hops(δL=>Matrix(complex(m)) for (δL, m)=ρ1) # convert to dense
    end

    H = Hamiltonian(Hops(), 0.0)

    function updateH!(H::Hamiltonian, ρ::AnyHops)
        verbose ? @info("Updating chemical potential for given filling.") : nothing
        H.h = ℋ_op(ρ) # get updated Hamiltonian
        H.μ = chemicalpotential(H.h, ks, filling; T=T)
        H
    end

    function update!(ρ1, ρ0)
        updateH!(H, ρ0)

        verbose ? @info("Updating the meanfield density matrix.") : nothing
        ϵ0 = densitymatrix!(ρ1, H.h, ks, H.μ; multimode=multimode, T=T, format=:dense) # get new meanfield and return the groundstate energy (density matrix was written to ρ1)

        if checkpoint != ""
            JLD.save(checkpoint, "mf", ρ1)
        end

        ϵ0 + ℋ_scalar(ρ1) # proper ground state energy
    end

    # Compute the ground state energy for the mean-field fixed point
    ϵ_GS, Error, converged = fixedpoint!(update!, ρ1, ρ0; iterations=iterations, tol=tol, kwargs...)

    if convergenceerror && !converged
        error("Convergence error.")
    end

    updateH!(H, ρ1)

    ρ1, ϵ_GS, H, converged, Error
end


# using JLD
import LatticeQM.Utils
import LatticeQM.Structure
import LatticeQM.Operators
import LatticeQM.Spectrum
import LatticeQM.TightBinding
using LatticeQM.Utils.Context

##############################################################################################################
##############################################################################################################
##############################################################################################################

import LatticeQM.TightBinding: Hops

#############################################
# Low-level TightBinding operator product
#############################################

using SparseArrays

function compute_AB_product(A_matrices::Dict{K,T1}, B_matrices::Dict{K,T2}) where {K,T1,T2}
    T3 = promote_type(T1, T2)
    AB_product = Dict{K,T3}()

    for R in keys(A_matrices)
        # AB_product_R = spzeros(eltype(T3), size(A_matrices[R])...)
        AB_product_R = zero(A_matrices[R])

        for R_prime in keys(B_matrices)
            R_diff = R .- R_prime  # Compute R - R'

            if haskey(A_matrices, R_diff)
                AB_product_R += A_matrices[R_diff] * B_matrices[R_prime]
            end
        end

        AB_product[R] = AB_product_R
    end

    return AB_product
end

SparseArrays.droptol!(H::Hops, args...) = SparseArrays.droptol!(H.data, args...)
function SparseArrays.droptol!(A::Dict{K,T}, tol::Real) where {K,T<:SparseMatrixCSC}
    for R in keys(A)
        SparseArrays.droptol!(A[R], tol)
    end
    return A
end


function compute_AB_product(A::Hops, B::Hops)
    Hops(compute_AB_product(A.data, B.data))
end

#############################################
# Low-level TightBinding operator product
#############################################

function initial_P_realspace(H, Emax, Emin, mu=0.0; beta=0.5)
    @assert Emax > mu > Emin "mu must be between Emin and Emax"
    alpha = min(beta / (Emax - mu), (1 - beta) / (mu - Emin))
    D = size(H, 1)

    id = Hops(Dict(TightBinding.zerokey(H) => complex(sparse(1.0 * I, D, D))))

    alpha * (mu * id - H) + beta * id
end

function canonicalpurification_grid_realspace(H, E_min, E_max, max_iter, filling=0.5; eps=1e-3, drop_tol=1e-4)

    N = size(H, 1)
    zk = TightBinding.zerokey(H)
    mu = real(tr(H[zk]) / N)

    P = initial_P_realspace(H, E_max, E_min, mu; beta=filling)
    P_2 = compute_AB_product(P, P)
    P_3 = compute_AB_product(P, P_2)

    

    converged = false
    prog = ProgressThresh(eps; desc="Density matrix search", showspeed=true)

    residual = Inf

    for iter in 1:max_iter
        c = abs(tr(P_2[zk] - P_3[zk]) / tr(P[zk] - P_2[zk]))

        @assert 0 <= c <= 1 "c must be between 0 and 1"

        if c >= 0.5
            P = ((1 + c) * P_2 - P_3) * (1 / c)
        else
            c < 0.5
            P = ((1 - 2 * c) * P + (1 + c) * P_2 - P_3) * (1 / (1 - c))
        end
        droptol!(P, drop_tol)

        P_2 = compute_AB_product(P, P)
        droptol!(P_2, drop_tol)
        P_3 = compute_AB_product(P, P_2)
        droptol!(P_2, drop_tol)

        residual = maximum(norm(P_2[R] - P[R], Inf) for R in keys(P))

        ProgressMeter.update!(prog, residual; showvalues=[(:iter, iter)])

        if residual < eps
            converged = true
            break # check for convergence
        end

    end
    ProgressMeter.finish!(prog)
    P_2 = nothing
    P_3 = nothing
    GC.gc()

    if !converged
        @warn "Maximum iterations reached without convergence or commutation (max residual = $residual)"
    end

    P
end

using Distributed 
import LatticeQM.Eigen 

function gridspectrumbound(H, ks)
    emin = min(real.(pmap(k -> Eigen.eigmin_sparse(H(k)), eachcol(ks)))...)
    emax = max(real.(pmap(k -> Eigen.eigmax_sparse(H(k)), eachcol(ks)))...)
    # emin = min(real.([LatticeQM.Eigen.eigmin_sparse(H(k)) for k in eachcol(ks)])...)
    # emax = max(real.([LatticeQM.Eigen.eigmax_sparse(H(k)) for k in eachcol(ks)])...)
    emin, emax
end

####### MWE
# lat2 = Geometries.honeycomb_twisted(6)
# H2 = Operators.graphene(lat2; tz=tz, format=:sparse, mode=:spinhalf)

# P = canonicalpurification_grid_realspace(H2, -3.9, 3.1, 20, 0.5; eps=1e-3)

##############################################################################################################
##############################################################################################################
##############################################################################################################  

# Specialized cases
solvehartreefock_purification(h::T, v, ρ_init, filling::Number, args...; kwargs...) where {T} = solveselfconsistent_purification!(ρ_init, HartreeFock(h, v), filling, args...; kwargs...)
solvefock_purification(h::T, v, ρ_init, filling::Number, args...; kwargs...) where {T} = solveselfconsistent_purification!(ρ_init, HartreeFock(h, v; hartree=false, fock=true), filling, args...; kwargs...)
solvehartree_purification(h::T, v, ρ_init, filling::Number, args...; kwargs...) where {T} = solveselfconsistent_purification!(ρ_init, HartreeFock(h, v; hartree=true, fock=false), filling, args...; kwargs...)

# Interface to solveselfconsistent!(ρ0, ...)
# solveselfconsistent_purification(ρ0, mf::MeanfieldGenerator, filling::Number, ks; kwargs...) = solveselfconsistent_purification!(ρ0, deepcopy(ρ0), mf, filling, ks; kwargs...)
# solveselfconsistent_purification(ρ0, mf::MeanfieldGenerator, filling::Number; klin, kwargs...) = solveselfconsistent_purification!(ρ0, deepcopy(ρ0), mf, filling; klin=klin, kwargs...)

# Interface to solveselfconsistent!(ρ0, ρ1, ...)
function solveselfconsistent_purification!(ρ0, mf::MeanfieldGenerator, filling::Float64, args...; kwargs...)
    sanitize!(ρ0)
    solveselfconsistent_purification!(ρ0, deepcopy(ρ0), mf, filling, args...; kwargs...)
end

# Interface to translate klin to solveselfconsistent!(ρ0, ρ1, ..., ks, ...)
solveselfconsistent_purification!(ρ0, ρ1, mf::MeanfieldGenerator, filling::Float64; klin::Int, kwargs...) = solveselfconsistent_purification!(ρ0, ρ1, mf, filling, Structure.regulargrid(nk=klin^2); kwargs...)

function solveselfconsistent_purification!(ρ0::T, ρ1::T, mf::MeanfieldGenerator, filling::Float64, ks; multimode=:global, kwargs...) where {T}
    c = trycontext(ensurecontext(multimode), SerialContext())
    solveselfconsistent_purification!(c, ρ0, ρ1, mf::MeanfieldGenerator, filling::Float64, ks; kwargs...)
end

using SharedArrays
import LatticeQM.Utils

solveselfconsistent_purification!(::SerialContext, ρ0::T, ρ1::T, args...; kwargs...) where {T} = solveselfconsistent_purification!(DummyContext(), ρ0, ρ1, args...; multimode=:serial, kwargs...)
function solveselfconsistent_purification!(::DistributedContext, ρ0::T, ρ1::T, args...; kwargs...) where {T}
    solveselfconsistent_purification!(DummyContext(), ρ0, ρ1, args...; multimode=:distributed, kwargs...)
end

function solveselfconsistent_purification!(::MultiThreadedContext, args...; kwargs...)
    error("Multithreaded context not implemented yet.")
end

"""
    solveselfconsistent!(ρ0, ρ1, ℋ_op, ℋ_scalar, filling, ks; convergenceerror=false, multimode=:serial, checkpoint::String="", hotstart=true, iterations=500, tol=1e-7, T=0.0, format=:dense, verbose::Bool=false, kwargs...)
    solveselfconsistent!(ρ0, ℋ_op, ℋ_scalar, filling, ks; kwargs...)
    solveselfconsistent!(ρ0, ρ1, ℋ_op, ℋ_scalar, filling; klin, kwargs...)
    solveselfconsistent!(ρ0, ℋ_op, ℋ_scalar, filling; klin, kwargs...)

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
function solveselfconsistent_purification!(::DummyContext, ρ0::T1, ρ1::T1, hartreefock::MeanfieldGenerator, filling::Float64, ks;
    convergenceerror=false, multimode=:serial, checkpoint::String="", callback=(x -> nothing), hotstart=true, iterations=500, tol=1e-7, T=0.0, format=:dense, verbose::Bool=false, kwargs...) where {T1}

    # if checkpoint != "" && isfile(checkpoint) && hotstart
    #     println("Loading checkpoint file as initial guess: $checkpoint")
    #     ρ0 = JLD.load(checkpoint, "mf")
    # end

    sanitize!(ρ0)
    sanitize!(ρ1)
    @assert TightBinding.ishermitian(ρ0) "Initial guess for density matrix must be hermitian."

    function update!(ρ1, ρ0)
        hartreefock(ρ0) # update meanfield (h is updated in-place)

        h0 = hMF(hartreefock)
        Emin, Emax = gridspectrumbound(h0, ks)
        DeltaE = Emax - Emin
        Emin, Emax = (Emin - 0.05 * DeltaE, Emax + 0.05 * DeltaE)
        purification_steps = 15

        ρ1 = canonicalpurification_grid_realspace(h0, Emin, Emax, purification_steps, filling; eps=1e-3, drop_tol=1e-4)

        ϵ0 = sum(tr(h0[R]'*ρ1[R]) for R in keys(h0)) # kinetic energy part (self-energy contribution is added below)

        callback(ρ1)

        # if checkpoint != ""
        #     verbose ? @info("Saving intermediate mean field...") : nothing
        #     JLD.save(checkpoint, "mf", DenseHops(ρ1))
        # end

        ϵ0 + hartreefock.ϵMF # proper ground state energy
    end

    # Compute the ground state energy for the mean-field fixed point
    ϵ_GS, residual, converged = fixedpoint!(update!, ρ1, ρ0; iterations=iterations, tol=tol, verbose=verbose, kwargs...)

    if convergenceerror && !converged
        error("Convergence error.")
    end
    hartreefock(ρ1) # update meanfield (h is updated in-place)

    ρ1, ϵ_GS, hartreefock, converged, residual
end 
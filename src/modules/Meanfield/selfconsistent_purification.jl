import LatticeQM.Utils
import LatticeQM.Structure
import LatticeQM.Operators
import LatticeQM.Spectrum
import LatticeQM.TightBinding
import LatticeQM.Parallel

import LatticeQM.TightBinding: Hops

#############################################
# Low-level TightBinding operator product
#############################################

using SparseArrays

function compute_AB_product(A_matrices::Dict{K,T1}, B_matrices::Dict{K,T2}) where {K,T1,T2}
    T3 = promote_type(T1, T2)
    AB_product = Dict{K,T3}()

    for R in keys(A_matrices)
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
# Canonical purification (real space)
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
            break
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

import LatticeQM.Eigen
import LatticeQM.Parallel

# Find spectral bounds by visiting each k once and computing both extremes —
# the original did two separate `pmap`s, doubling the work. kspace_foreach!
# distributes the loop across the chosen executor (auto-picks distributed
# if workers exist, else threaded, else serial).
function gridspectrumbound(H, ks)
    exec = Parallel.to_executor(:auto)
    Parallel.configure_blas!(exec; verbose=false)

    emin = Ref(Inf)
    emax = Ref(-Inf)
    bounds_lock = Threads.SpinLock()

    Parallel.kspace_foreach!(ks, exec) do _scratch, _j, k
        Hk = H(k)
        a = real(Eigen.eigmin_sparse(Hk))
        b = real(Eigen.eigmax_sparse(Hk))
        Base.@lock bounds_lock begin
            emin[] = min(emin[], a)
            emax[] = max(emax[], b)
        end
    end
    emin[], emax[]
end

##############################################################################
# SCF interfaces
##############################################################################

# Specialized cases
solvehartreefock_purification(h::T0, v, ρ_init, filling::Number, args...; kwargs...) where {T0} = solveselfconsistent_purification!(ρ_init, HartreeFock(h, v), filling, args...; kwargs...)
solvefock_purification(h::T0, v, ρ_init, filling::Number, args...; kwargs...) where {T0} = solveselfconsistent_purification!(ρ_init, HartreeFock(h, v; hartree=false, fock=true), filling, args...; kwargs...)
solvehartree_purification(h::T0, v, ρ_init, filling::Number, args...; kwargs...) where {T0} = solveselfconsistent_purification!(ρ_init, HartreeFock(h, v; hartree=true, fock=false), filling, args...; kwargs...)

function solveselfconsistent_purification!(ρ0, mf::MeanfieldGenerator, filling::Float64, args...; kwargs...)
    sanitize!(ρ0)
    solveselfconsistent_purification!(ρ0, deepcopy(ρ0), mf, filling, args...; kwargs...)
end

solveselfconsistent_purification!(ρ0, ρ1, mf::MeanfieldGenerator, filling::Float64; klin::Int, kwargs...) = solveselfconsistent_purification!(ρ0, ρ1, mf, filling, Structure.regulargrid(nk=klin^2); kwargs...)

function solveselfconsistent_purification!(ρ0::T0, ρ1::T0, mf::MeanfieldGenerator, filling::Float64, ks; multimode=:auto, kwargs...) where {T0}
    mm = (multimode === :global) ? :auto : multimode
    exec = Parallel.to_executor(mm)

    if exec isa Parallel.DistributedExec
        return _solveselfconsistent_purification_impl!(ρ0, ρ1, mf, filling, ks; multimode=:distributed, kwargs...)
    elseif exec isa Parallel.ThreadedExec
        # Same story as solveselfconsistent!: threading was a stub before
        # Parallel; now the downstream `Parallel.kspace_*` accumulators handle
        # it correctly.
        return _solveselfconsistent_purification_impl!(ρ0, ρ1, mf, filling, ks; multimode=:multithreaded, kwargs...)
    else
        return _solveselfconsistent_purification_impl!(ρ0, ρ1, mf, filling, ks; multimode=:serial, kwargs...)
    end
end

"""
    solveselfconsistent_purification!(ρ0, ρ1, mf, filling, ks; kwargs...)

Real-space canonical-purification analogue of [`solveselfconsistent!`](@ref).
Each SCF step replaces k-space diagonalization with a fixed-step McWeeny-style
purification of the current mean-field Hamiltonian. Useful for very large
systems where per-k diagonalization is the bottleneck.

Returns `(ρ, ϵ_GS, mf, converged, residual)`. Same kwargs as
`solveselfconsistent!`.
"""
function _solveselfconsistent_purification_impl!(ρ0::T1, ρ1::T1, hartreefock::MeanfieldGenerator, filling::Float64, ks;
    multimode=:serial,
    purification_steps::Int=15,
    purification_eps::Real=1e-3,
    purification_drop_tol::Real=1e-4,
    kwargs...) where {T1}

    update_ρ!(ρ1, hf) = begin
        h0 = hMF(hf)
        Emin, Emax = gridspectrumbound(h0, ks)
        ΔE = Emax - Emin
        Emin, Emax = (Emin - 0.05 * ΔE, Emax + 0.05 * ΔE)

        ρ_new = canonicalpurification_grid_realspace(h0, Emin, Emax, purification_steps, filling;
                                                     eps=purification_eps, drop_tol=purification_drop_tol)

        # Write back into the caller's ρ1 in place. Re-binding `ρ1 = ρ_new`
        # would only update the local parameter — `fixedpoint!`'s outer view
        # would never see the new density and the iteration would converge
        # immediately on the initial guess.
        for k in keys(ρ_new)
            ρ1[k] = ρ_new[k]
        end

        sum(tr(h0[R]' * ρ1[R]) for R in keys(h0))
    end

    _scf_driver!(ρ0, ρ1, hartreefock, update_ρ!; kwargs...)
end

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

# Real-space convolution restricted to keys(A): (A * B)(R) = Σ_{R'} A(R-R') B(R')
# for every R already in keys(A). This is intentional — the McWeeny iteration
# would otherwise grow the support exponentially per step, which is fatal on
# small unit cells where droptol! can't prune the dense O(D²) blocks fast
# enough. Practical implication: P stays pinned to H's sparsity pattern, so
# the converged ρ is the projection of the true density matrix onto that
# pattern. For accurate ρ, ensure the H you pass already covers the support
# you need — e.g. by precomputing `H * H' * ...` enough times to grow keys
# before invoking purification.
function compute_AB_product(A_matrices::Dict{K,T1}, B_matrices::Dict{K,T2}) where {K,T1,T2}
    T3 = promote_type(T1, T2)
    AB_product = Dict{K,T3}()

    for R in keys(A_matrices)
        AB_product_R = zero(A_matrices[R])

        for R_prime in keys(B_matrices)
            R_diff = R .- R_prime
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
    # Drop small entries within each block, then prune blocks that became
    # all-zero so the convolution support doesn't grow without bound.
    to_remove = K[]
    for R in keys(A)
        SparseArrays.droptol!(A[R], tol)
        nnz(A[R]) == 0 && push!(to_remove, R)
    end
    for R in to_remove
        delete!(A, R)
    end
    return A
end


function compute_AB_product(A::Hops, B::Hops)
    Hops(compute_AB_product(A.data, B.data))
end

# ‖A[R] - B[R]‖_∞, treating missing blocks as zero. Used by McWeeny
# convergence check after the convolution fix made A and B's key sets diverge.
function _block_diff_norm(A::Hops, B::Hops, R)
    if haskey(A, R) && haskey(B, R)
        return norm(A[R] - B[R], Inf)
    elseif haskey(A, R)
        return norm(A[R], Inf)
    elseif haskey(B, R)
        return norm(B[R], Inf)
    else
        return 0.0
    end
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
        droptol!(P_3, drop_tol)

        # Convergence: P should equal P² (idempotent projector). Use the
        # union of supports because droptol! can prune blocks asymmetrically
        # between P and P_2 (e.g. one becomes all-zero and is removed).
        residual = maximum(_block_diff_norm(P, P_2, R) for R in union(keys(P), keys(P_2)))

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

Real-space canonical-purification analogue of [`solveselfconsistent`](@ref).
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
    residual_norm::Symbol=:density,
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

    # See `_solveselfconsistent_impl!` for the rationale on residual choice.
    # Purification grows `keys(ρ)` over iterations, so commutator residual
    # *can* converge here (unlike the k-space SCF with truncated `ρ`), but
    # only if McWeeny is run to a tight enough tolerance per step.
    residual_fn = if residual_norm === :density
        nothing
    elseif residual_norm === :commutator
        (ρ1_, ρ0_) -> commutator_kspace_norm(hMF(hartreefock), ρ0_, ks;
                                              multimode=multimode)
    else
        throw(ArgumentError("residual_norm must be :density or :commutator, got $(residual_norm)"))
    end

    _scf_driver!(ρ0, ρ1, hartreefock, update_ρ!;
                 residual_fn=residual_fn, kwargs...)
end

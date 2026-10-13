using ProgressMeter

"""
    fixedpoint!(f!, x1, x0; iterations=500, tol=1e-7, β=1.0, p_norm::Real=2,
                relative=true, acceleration=:linear, anderson_depth=5,
                show_trace=false, clear_trace=false, verbose=true,
                log_callback=nothing, residual_fn=nothing)

Damped fixed-point iteration with optional Anderson (Type-II) acceleration.
`f!(x1, x0)` should overwrite `x1` with the next iterate computed from `x0`
and may return a scalar (e.g. ground-state energy at this step).

After each call to `f!`:

  * **Linear mixing** (`acceleration=:linear`, default): the new iterate is
    `β * f(x_k) + (1 - β) * x_k`. `β=1` is undamped.
  * **Anderson acceleration** (`acceleration=:anderson`, opt-in): keeps the
    last `anderson_depth` `(x_i, f(x_i))` pairs and solves a small
    least-squares problem to combine them. `β` acts as a damping parameter
    (`β=1` is undamped Anderson). Falls back to plain mixing on the first
    step (no history yet) and whenever the LS solve produces a non-finite
    update. Typical SCF problems converge ~5–10× faster than linear mixing
    on well-conditioned interacting systems. **Caveat:** at `U ≈ 0` and
    `T = 0` the residuals between iterations are nearly collinear, the LS
    becomes ill-conditioned, and γ blows up — Anderson is opt-in rather
    than the default for that reason. A future version may add automatic
    rank-revealing QR pruning to make it safe by default.

The residual is the raw fixed-point defect `‖f(x0) - x0‖`, divided by
`‖f(x0)‖` when `relative=true` and the norm is nonzero. It is measured
before damping or Anderson acceleration, so a small step cannot hide an
unconverged fixed point. An optional `residual_fn=(x1, x0) -> Real` adds a
second convergence condition; the reported residual is the larger of the
two. The callback sees the raw `(f(x0), x0)` pair before mixing.
Once converged, both buffers hold the validated `x0`; no further mixing or
acceleration is applied to it.
Set `log_callback=(iter, ϵ, residual) -> ...` to log every step.

Sanity-checked against `f_a(x) = 1/2 * (a/x + x)`, fixed point `√a`.
"""
function fixedpoint!(f!, x1, x0;
    iterations=500,
    tol=1e-7, β=1.0,
    p_norm=2,
    relative::Bool=true,
    acceleration::Symbol=:linear,
    anderson_depth::Int=5,
    show_trace=false,
    clear_trace=false,
    verbose=true,
    log_callback=nothing,
    residual_fn=nothing,
    )

    use_anderson = acceleration === :anderson && anderson_depth > 0
    use_anderson || acceleration === :linear ||
        throw(ArgumentError("acceleration must be :linear or :anderson, got $(acceleration)"))

    converged = false
    ϵ0 = 0.0

    if show_trace
        prog = ProgressThresh(tol; desc="FIXPOINT SEARCH ", showspeed=true)
    end

    # Anderson history setup: stable key ordering + flattened scratch buffers.
    if use_anderson
        key_order = sort!(collect(keys(x0)))
        key_sizes = Int[length(x0[k]) for k in key_order]
        n_var = sum(key_sizes)
        T_el = eltype(first(values(x0)))
        flat_x = Vector{T_el}(undef, n_var)
        flat_f = Vector{T_el}(undef, n_var)
        hist_x = Vector{Vector{T_el}}()    # x_{k-m}, ..., x_k
        hist_f = Vector{Vector{T_el}}()    # f(x_{k-m}), ..., f(x_k)
    end

    residual = Inf
    iter = 0
    while iter < iterations
        iter += 1

        if use_anderson
            _flatten!(flat_x, x0, key_order, key_sizes)
        end

        ϵ0 = f!(x1, x0)

        # Pair blocks by key, including keys introduced by the update.
        _sync_keysets!(x1, x0)
        diff_norm = norm([norm(x1[L] .- x0[L], p_norm) for L in keys(x0)], p_norm)
        x1_norm = relative ? norm([norm(x1[L], p_norm) for L in keys(x0)], p_norm) : 0.0
        residual = relative && x1_norm > eps() ? diff_norm / x1_norm : diff_norm
        if residual_fn !== nothing
            residual = max(residual, residual_fn(x1, x0))
        end

        if residual < tol
            # The measured defect belongs to x0. An Anderson extrapolation
            # can leave this fixed point even when the defect is already small.
            converged = true
            for L in keys(x0)
                x1[L] .= x0[L]
            end
        elseif use_anderson
            if sort!(collect(keys(x0))) != key_order
                # Layout changed: stale flat history has the wrong length and
                # would mispair blocks. Rebuild buffers and restart the
                # acceleration history from this iterate.
                key_order = sort!(collect(keys(x0)))
                key_sizes = Int[length(x0[k]) for k in key_order]
                n_var = sum(key_sizes)
                flat_x = Vector{T_el}(undef, n_var)
                flat_f = Vector{T_el}(undef, n_var)
                empty!(hist_x); empty!(hist_f)
                _flatten!(flat_x, x0, key_order, key_sizes)
            end
            _flatten!(flat_f, x1, key_order, key_sizes)

            push!(hist_x, copy(flat_x))
            push!(hist_f, copy(flat_f))
            while length(hist_x) > anderson_depth + 1
                popfirst!(hist_x); popfirst!(hist_f)
            end

            flat_next = _anderson_update(flat_x, flat_f, hist_x, hist_f, β)

            _unflatten!(x1, flat_next, key_order, key_sizes)
            # The unconstrained LS doesn't know about the Hops Hermiticity
            # constraint (x1[L]' == x1[-L]); project back. Cheap, and
            # required for downstream meanfield assertions.
            _hermitianize_if_supported!(x1)
            for k in key_order
                x0[k] .= x1[k]
            end

        else
            # Linear damped mixing (key sets are synced above, so same-key
            # pairing over keys(x0) covers everything)
            for δL in keys(x0)
                @. x1[δL] = β * x1[δL] + (1 - β) * x0[δL]
            end

            for δL in keys(x0)
                @. x0[δL] = x1[δL]
            end
        end

        if show_trace
            ProgressMeter.update!(prog, residual)
        end

        if log_callback !== nothing
            log_callback(iter, ϵ0, residual)
        end

        converged && break
    end

    if show_trace && clear_trace
        ProgressMeter.finish!(prog)
    end

    if verbose
        @info("Total #iterations: $iter   Residual error: $residual")
    end

    if !converged
        @warn("Did not converge to requested target tolerance. Residual: $residual")
    end

    ϵ0, residual, converged
end

# ----------------------------------------------------------------------------
# Anderson Type-II update with damping β. Returns the next iterate as a flat
# Vector. The first step (no usable history) and any LS solve that produces a
# non-finite γ both fall back to plain damped mixing.
# ----------------------------------------------------------------------------
function _anderson_update(flat_x::AbstractVector{T}, flat_f::AbstractVector{T},
                          hist_x::Vector{Vector{T}}, hist_f::Vector{Vector{T}},
                          β::Real) where {T}
    m = length(hist_x) - 1
    if m == 0
        return @. β * flat_f + (1 - β) * flat_x
    end

    n = length(flat_x)
    ΔX = Matrix{T}(undef, n, m)
    ΔF = Matrix{T}(undef, n, m)
    @inbounds for i in 1:m
        @views ΔX[:, i] .= hist_x[i+1] .- hist_x[i]
        @views ΔF[:, i] .= hist_f[i+1] .- hist_f[i]
    end
    ΔG = ΔF .- ΔX                # residual differences
    g_k = flat_f .- flat_x       # current residual

    γ = try
        ΔG \ g_k                 # least-squares solve, m unknowns
    catch
        return @. β * flat_f + (1 - β) * flat_x
    end

    if !all(isfinite, γ)
        return @. β * flat_f + (1 - β) * flat_x
    end

    # x_{k+1} = (1-β)(x_k - ΔX γ) + β (f_k - ΔF γ)
    flat_next = similar(flat_x)
    @views @. flat_next = (1 - β) * flat_x + β * flat_f
    BLAS_axpy_combination!(flat_next, ΔX, ΔF, γ, β)
    flat_next
end

# In-place: flat_next .-= ((1-β)*ΔX + β*ΔF) * γ
function BLAS_axpy_combination!(flat_next::AbstractVector, ΔX::AbstractMatrix,
                                ΔF::AbstractMatrix, γ::AbstractVector, β::Real)
    @inbounds for i in eachindex(flat_next)
        s = zero(eltype(flat_next))
        @inbounds for j in eachindex(γ)
            s += ((1 - β) * ΔX[i, j] + β * ΔF[i, j]) * γ[j]
        end
        flat_next[i] -= s
    end
    flat_next
end

# ----------------------------------------------------------------------------
# Flatten / unflatten helpers for Hops-like containers (Dict-of-matrix).
# ----------------------------------------------------------------------------

# Extend both iterates to the union of their block key sets, filling missing
# blocks with zeros of the partner's shape/backing. Returns true if either
# container gained a key. Keeps mixing/residual/copy-back well-defined when
# an update adds blocks.
function _sync_keysets!(x1, x0)
    changed = false
    for δL in keys(x0)
        haskey(x1, δL) && continue
        x1[δL] = 0 .* x0[δL]
        changed = true
    end
    for δL in keys(x1)
        haskey(x0, δL) && continue
        x0[δL] = 0 .* x1[δL]
        changed = true
    end
    changed
end

function _flatten!(out::AbstractVector, x, key_order, key_sizes)
    offset = 0
    @inbounds for (k, sz) in zip(key_order, key_sizes)
        block = x[k]
        @views out[offset+1:offset+sz] .= vec(block)
        offset += sz
    end
    out
end

_hermitianize_if_supported!(x) = x  # generic fallback
function _hermitianize_if_supported!(x::TightBinding.Hops)
    TightBinding.hermitianize!(x)
    x
end

function _unflatten!(x, v::AbstractVector, key_order, key_sizes)
    offset = 0
    @inbounds for (k, sz) in zip(key_order, key_sizes)
        block = x[k]
        block .= reshape(@view(v[offset+1:offset+sz]), size(block))
        offset += sz
    end
    x
end

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

The default residual is `‖x1 - x0‖ / max(‖x1‖, eps())` if `relative=true`,
else absolute. Pass `residual_fn=(x1, x0) -> Real` to override it — the
callback is evaluated immediately after `f!`, before mixing rewrites `x0`,
so it sees the raw `(x_{k+1}^{f}, x_k)` pair. The Anderson least-squares
math always uses the iterate-difference form internally regardless of
`residual_fn`. Set `log_callback=(iter, ϵ, residual) -> ...` to log every
step without enabling the progress bar.

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

        # Custom residual must be evaluated before mixing overwrites x0 (and
        # before the Anderson combination overwrites x1). It typically
        # depends on auxiliary state populated by `f!` — e.g. the SCF driver
        # caches `H_MF[ρ_0]` on the meanfield generator and computes the
        # commutator `‖[H_MF, ρ_0]‖` here.
        custom_residual = residual_fn === nothing ? nothing : residual_fn(x1, x0)

        if use_anderson
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
            # Refresh flat_next from the projected x1 so the residual norm
            # below reflects what the next iteration actually sees.
            _flatten!(flat_next, x1, key_order, key_sizes)
            for k in key_order
                x0[k] .= x1[k]
            end

            diff_norm = norm(flat_next .- flat_x, p_norm)
            if relative
                x1_norm = norm(flat_next, p_norm)
                residual = x1_norm > eps() ? diff_norm / x1_norm : diff_norm
            else
                residual = diff_norm
            end
        else
            # Linear damped mixing
            for δL in keys(x0)
                @. x1[δL] = β * x1[δL] + (1 - β) * x0[δL]
            end

            diff_norm = norm(values(x1) .- values(x0), p_norm)
            if relative
                x1_norm = norm(values(x1), p_norm)
                residual = x1_norm > eps() ? diff_norm / x1_norm : diff_norm
            else
                residual = diff_norm
            end

            for δL in keys(x0)
                @. x0[δL] = x1[δL]
            end
        end

        if custom_residual !== nothing
            residual = custom_residual
        end

        if show_trace
            ProgressMeter.update!(prog, residual)
        end

        if log_callback !== nothing
            log_callback(iter, ϵ0, residual)
        end

        if residual < tol
            converged = true
            break
        end
    end

    if show_trace && clear_trace
        ProgressMeter.finish!(prog)
    end

    if verbose
        @info("Total #iterations: $iter   Residual error: $residual")
    end

    if !converged
        @warn("Did not convergence to requested target tolerance. Residual: $residual")
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

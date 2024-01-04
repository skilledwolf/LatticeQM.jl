# using Printf

using ProgressMeter

# abstract type AbstractFixedPointState{T} end

# struct FixedPointState{T} <: AbstractFixedPointState{T}
#     x_current::T
#     x_proposed::T
#     f_current::Float64
#     f_proposed::Float64
#     residual::Float64
#     converged::Bool
# end

# @with_kw struct FixedPointSettings
#     max_iterations::Int = 500 
#     tol::Float64 = sqrt(eps())
#     β::Float64 = 1.0
#     p_norm::Real = 2
#     show_trace::Bool = false
#     verbose::Bool = false
# end

# function update!(state::AbstractFixedPointState, settings)
#     error("Should be implemented by the user.")
# end

# function step!(state::FixedPointState{<:AbstractMatrix}, settings)
#     @. state.x_proposed = settings.β * state.x_proposed + (1 - settings.β) * state.x_current
#     state.residual = norm(state.x_proposed .- state.x_current, settings.p_norm) / norm(state.x_current, settings.p_norm)
#     state.x_current .= state.x_proposed
# end

# function iterate!(state::AbstractFixedPointState, settings::FixedPointSettings)
#     state.converged = false

#     settings.show_trace &&  prog = ProgressThresh(settings.tol, "FIXPOINT ITERATION "; showspeed=true)

#     iter = 0
#     while iter < iteration.max_iterations && state.residual > iteration.tol
#         iter += 1
#         update!(state, settings)
#         step!(state, settings)

#         settings.show_trace &&  ProgressMeter.update!(prog, state.residual)
#     end

#     state.converged = state.residual <= settings.tol

#     verbose && @info("Total #iterations: $iter   Residual error: $(state.residual)")
#     !state.converged && @warn("Did not convergence to requested target tolerance. Residual: $(state.residual)")

#     state
# end


"""
    fixedpoint!(f!, x1, x0; iterations=500, tol=1e-7, β=1.0, p_norm::Real=2, show_trace=false, clear_trace=false)

Performs a simple fixed point iteration. The function f!(x1,x0) should override with x1 with x2 and x0 with x1 and
can optionally return a scalar value (for example ground state energy at iteration step). 

Fixedpoint iteration, tested on the square-root example
f_a(x) = 1/2 * (a/x+x)
which has the fixed point x0 = sqrt(a).
"""
function fixedpoint!(f!, x1, x0;
    iterations=500,
    tol=1e-7, β=1.0,
    p_norm=2,#Inf,
    show_trace=false,
    verbose=true
    )

    converged = false
    ϵ0 = 0.0

    if show_trace
        prog = ProgressThresh(tol; desc="FIXPOINT SEARCH ", showspeed=true)
    end

    residual = 1.0
    iter = 0
    while iter < iterations
        iter += 1

        # Perform a step
#         timediff = @timed ϵ0 = f!(x1, x0)
        ϵ0 = f!(x1, x0)

        # Convergence acceleration ("damped fixed point iterater")
        # The new x0 for the next step is:
        for δL = keys(x0)
            @. x1[δL] .= β * x1[δL] + (1 - β) * x0[δL]
        end

        # Convergence?
        # residual = norm(values(x1) .- values(x0), p_norm)
        residual = norm(values(x1) .- values(x0), p_norm)/norm(values(x1), p_norm)

        # Convergence acceleration ("damped fixed point iterater")
        # The new x0 for the next step is:
        for δL = keys(x0)
            @. x0[δL] .= x1[δL]
        end

        if show_trace
            ProgressMeter.update!(prog, residual)
        end

        if residual < tol
            converged = true
            break
        end

        if iter == iterations
            break
        end
    end

    if verbose
        @info("Total #iterations: $iter   Residual error: $residual")
    end

    if !converged
        @warn("Did not convergence to requested target tolerance. Residual: $residual")
    end

    ϵ0, residual, converged
end



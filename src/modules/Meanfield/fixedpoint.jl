using Printf

using ProgressMeter

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
    p_norm=Inf,
    show_trace=false,
    clear_trace=false,
    verbose=true
    )

    converged = false
    ϵ0 = 0.0

    if show_trace
        # println("New step.")
        prog = ProgressThresh(tol, "FIXPOINT SEARCH ")
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
            ProgressMeter.update!(prog, error, showvalues = [(:iter,iter), (:time, t1)])
        end

        if error < tol
            converged = true
            break
        end

        if iter == iterations
            break
        end

        # Convergence acceleration ("damped fixed point iterater")
        # The new x0 for the next step is:
        for δL=keys(x0)
            @. x0[δL][:] .= β * x1[δL][:] + (1-β) * x0[δL][:]
        end
    end

    if verbose
        @info("Total #iterations: $iter   Residual error: $error")
    end

    if !converged
        @warn("Did not convergence to requested target tolerance. Residual: $error")
    end

    ϵ0, error, converged
end



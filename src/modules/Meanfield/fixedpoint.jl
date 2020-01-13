using Printf

function fixedpoint!(f!, x1, x0;
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


###################################################################################################
# Backwards compatibility
###################################################################################################
@legacyalias fixedpoint search_fixedpoint


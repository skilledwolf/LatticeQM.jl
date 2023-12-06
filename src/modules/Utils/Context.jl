module Context

import Distributed#, Threads

export ExecutionContext, SerialContext, MultiThreadedContext, DistributedContext, DummyContext
export setglobalcontext, isavailable, trycontext, ensurecontext

abstract type ExecutionContext end

struct SerialContext <: ExecutionContext end
struct DummyContext <: ExecutionContext end
struct MultiThreadedContext <: ExecutionContext end
struct DistributedContext <: ExecutionContext end
# struct AutomaticContext <: ExecutionContext end

ensurecontext(context::ExecutionContext) = context
function ensurecontext(context::Symbol=:auto)
    if context == :serial
        return SerialContext()
    elseif context == :distributed 
        return DistributedContext()
    elseif context == :multithread
        return MultiThreadedContext()
    elseif context == :auto
        return getautocontext()
    elseif context == :global
        return getglobalcontext()
    else
        @warn "Unknown context keyword. Defaulting to :serial."
        return SerialContext()
    end
end

isavailable(::SerialContext) = true
# isavailable(::AutomaticContext) = true
isavailable(context::DistributedContext) = (Distributed.nworkers() > 1)
isavailable(context::MultiThreadedContext) = (Threads.nthreads() > 1)

function getautocontext()
    if isavailable((distr = DistributedContext(); distr))
        return distr
    elseif isavailable((multith = MultiThreadedContext(); multith))
        return multith
    else
        return SerialContext()
    end
end

# Global variable for setting user preference
global user_preferred_context = getautocontext()::ExecutionContext

# Function to set global context
setglobalcontext(context::Symbol=:auto) = (global user_preferred_context = ensurecontext(context))
setglobalcontext(context::ExecutionContext) = (global user_preferred_context = context)
getglobalcontext() = user_preferred_context

function trycontext(context::ExecutionContext, fallbackcontext::ExecutionContext=SerialContext(); raise=false, warn=true)
    avail = isavailable(context) 

    if avail 
        return context
    else
        if raise 
            error("Execution context is not available: ", string(typeof(context)))
        end
        if warn
            @warn "The requested context is not available: $(string(typeof(context)))"
        end
        return fallbackcontext
    end
end

end
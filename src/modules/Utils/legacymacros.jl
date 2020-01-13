macro legacyalias(f, oldname)
"""
This is macro is takes function f(...) and provides an alias oldname(...) with a
deprecation warning.

Usage:
newfunc(x) = x

@legacyalias(newfunc, oldfunc)
"""
    return quote
        name1 = $(string(oldname))
        name2 = $(string(f))
        function $(esc(oldname))(args...; kwargs...)
            @warn("$name1(...) was renamed to $name2(...) and is marked for removal.")
            $(esc(f))(args...; kwargs...)
        end

    end
end

macro legacymoved(f, newname)
"""
This is macro raises an error when function f(...) is used and states the new location of
the definition.

@legacymoved(oldfunc, "Newmodule.newfunc")
"""
    return quote
        name1 = $(string(f))
        name2 = $(string(newname))
        function $(esc(f))(args...; kwargs...)
            error("$name1(...) was moved to $name2(...).")
        end

    end
end

macro legacyremoved(f)
"""
This is macro raises an error when function f(...) is used.

@legacyremoved(oldfunc)
"""
    return quote
        name1 = $(string(f))
        function $(esc(f))(args...; kwargs...)
            error("$name1(...) was removed.")
        end

    end
end
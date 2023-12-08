

function getoperator(lat::Lattice, name::String, args...; kwargs...)

    if name == "MX" || name == "SX" || name == "sx"
        return getsx(lat, args...; kwargs...)
    elseif name == "MY" || name == "SY" || name == "sy"
        return getsy(lat, args...; kwargs...)
    elseif name == "MZ" || name == "SZ" || name == "sz" || name == "spin"
        return getsz(lat, args...; kwargs...)
    elseif name == "spinUP" || name == "spinup"
        return getsup(lat, args...; kwargs...)
    elseif name == "spinDOWN"  || name == "spindown"
        return getsdown(lat, args...;kwargs...)
    elseif name =="Sn"
        return Sn(lat, args...; kwargs...)
    elseif name =="layer"
        return layerprojection(lat, args...; kwargs...)
    elseif name == "sublattice"
        return sublatticeprojection(lat, args...; kwargs...)
    elseif name == "sublatticeA"
        return sublatticeprojection(lat, "A", args...; kwargs...)
    elseif name == "sublatticeB"
        return sublatticeprojection(lat, "B", args...; kwargs...)
    elseif name == "sublatticeAspin"
        return sublatticeprojection(lat, "A", 2, args...; kwargs...)
    elseif name == "sublatticeBspin"
        return sublatticeprojection(lat, "B", 2, args...; kwargs...)
    elseif name == "Hubbard"
        error("Load operator '$name' from Operators.gethubbard(...).")
    elseif name == "CappedYukawa"
        error("Load operator '$name' from Meanfield.getcappedyukawa(...).")
    else
        error("Requested operator '$name' not found.")
    end
end

getoperator(lat::Lattice, names::AbstractVector{String}, args...; kwargs...) = [getoperator(lat, name, args...; kwargs...) for name=names]

# import ..TightBinding: expvalf
# getprojector(lat::Lattice, name::String, args...; kwargs...) = expvalf(getoperator(lat, name, args...; kwargs...))
# getprojector(lat::Lattice, names::AbstractVector{String}, args...; kwargs...) = [expvalf(getoperator(lat, name, args...; kwargs...)) for name=names]

function getoperator(lat::Lattice, name::String, args...; kwargs...)

    if name == "MX" || name == "SX"
        return SX(lat, args...; kwargs...)
    elseif name == "MY" || name == "SY"
        return SY(lat, args...; kwargs...)
    elseif name == "MZ" || name == "SZ" || name == "spin"
        return SZ(lat, args...; kwargs...)
    elseif name == "spinUP"
        return Sup(lat, args...; kwargs...)
    elseif name == "spinDOWN"
        return Sdown(lat, args...;kwargs...)
    elseif name =="Sn"
        return S_n(lat, args...; kwargs...)
    elseif name =="layer"
        return layer_op(lat, args...; kwargs...)
    elseif name == "sublattice"
        return sublattice_N(lat, args...; kwargs...)
    elseif name == "sublatticeA"
        return sublattice_N(lat, 0, args...; kwargs...)
    elseif name == "sublatticeB"
        return sublattice_N(lat, 1, args...; kwargs...)
    elseif name == "sublatticeAspin"
        return sublattice_N(lat, 0, 2, args...; kwargs...)
    elseif name == "sublatticeBspin"
        return sublattice_N(lat, 1, 2, args...; kwargs...)
    elseif name == "Hubbard"
        error("Load operator '$name' from Meanfield.gethubbard(...).")
    elseif name == "CappedYukawa"
        error("Load operator '$name' from Meanfield.getcappedyukawa(...).")
    else
        error("Requested operator '$name' not found.")
    end
end

getoperator(lat::Lattice, names::AbstractVector{String}, args...; kwargs...) = [getoperator(lat, name, args...; kwargs...) for name=names]

getprojector(lat::Lattice, name::String, args...; kwargs...) = expvalf(getoperator(lat, name, args...; kwargs...))
getprojector(lat::Lattice, names::AbstractVector{String}, args...; kwargs...) = [expvalf(getoperator(lat, name, args...; kwargs...)) for name=names]

# Note: getprojector is a "bad" name and should be moved to get_expvalf
#       generally I use mostly for expectation values in bandstructures,
#       not projections onto basis

################################################################################
################################################################################
################################################################################

function expvalf(𝑶::AbstractMatrix)


    f(k, ψ, ϵ) = real.(ψ' * 𝑶 * ψ)

    f
end

function expvalf(𝑶::Function)

    f(k, ψ, ϵ) = real.(ψ' * 𝑶(k) * ψ)

    f
end

################################################################################
################################################################################
################################################################################


latticemap(f,lat) = map(f, allpositions(lat))

function layer_op(lat::Lattice, d=1)
    """
    typically n=0 for sublattice A and n=1 for sublattice B and d=2 for spin-1/2
    """

    z = extrapositions(lat, "z")

    zmin = minimum(z)
    zmax = maximum(z)

    z = 2.0 .* ( (z .- zmin) ./ (zmax-zmin) .- 0.5 )

    kron(Diagonal(z), Diagonal(ones(d)))
end

function sublattice_N(lat::Lattice, n, d=1)
    """
    typically n=0 for sublattice A and n=1 for sublattice B and d=2 for spin-1/2
    """

    sublattice = extrapositions(lat, "sublattice")

    filtered_sublattice = Array{Float64}(sublattice .== n)


    kron(Diagonal(filtered_sublattice), Diagonal(ones(d)))
end


###################################################################################################
# Backwards compatibility
###################################################################################################
export get_operator
@legacyalias getoperator get_operator

export get_projector
@legacyalias getprojector get_projector

@legacymoved initial_guess "Meanfield.initialguess"

export set_filling!
@legacymoved set_filling! "Materials.setfilling!"

export add_chemicalpotential!
@legacymoved add_chemicalpotential! "Materials.addchemicalpotential!"
getzeeman(args...; kwargs...) = addzeeman!(Hops(), args...; kwargs...)

function addzeeman!(hops, lat::Lattice, Mv::Function)
    zero0 = zeros(Int, latticedim(lat))

    N = countorbitals(lat)
    R = allpositions(lat)

    σn(vec) = sum(vec[i] .* σs[i] for i=1:3)

    mat = spzeros(ComplexF64, 2*N, 2*N)
    for i_=1:N
        mat[1+2*(i_-1):2+2*(i_-1),1+2*(i_-1):2+2*(i_-1)] .= σn(Mv(R[:,i_]))
    end

    addhops!(hops, Hops(zero0 => mat))
end

function addzeeman!(hops, lat::Lattice, Mv::Vector{Float64}; format=:dense)
    # Only go through the trouble of constructing this matrix for finite Mv
    if isapprox(norm(Mv), 0; atol=sqrt(eps()))
        return 0.0
    end

    zero0 = zeros(Int, latticedim(lat))

    N = countorbitals(lat)
    R = positions(lat)

    σn = sum(Mv[i] .* σs[i] for i=1:3)

    newhops = Dict( zero0 => kron(Matrix(1.0I, N, N), σn) )
    addhops!(hops, newhops)

    hops
end
addzeeman!(hops, lat::Lattice, M0::Float64; kwargs...) = addzeeman!(hops, lat::Lattice, [0.0,0.0,M0]; kwargs...)



###################################################################################################
# Backwards compatibility
###################################################################################################
export zeeman_hops
@legacyalias getzeeman zeeman_hops

export add_zeeman!
@legacyalias addzeeman! add_zeeman!

const LatticeHopsAbstract = Dict{Vector{Int}, <:AbstractMatrix{ComplexF64}}
const LatticeHopsDense = Dict{Vector{Int}, <:Matrix{ComplexF64}}
const LatticeHopsSparse = Dict{Vector{Int}, <:SparseMatrixCSC{ComplexF64}}
const LatticeHops = Union{LatticeHopsDense, LatticeHopsSparse, LatticeHopsAbstract}

# LatticeHopsSparse{T}(args...; kwargs...) where T<:SparseMatrixCSC{<:Complex} = Dict{Vector{Int}, T}(Vararg{Pair,N} where N)#Dict{Vector{Int}, T}(args...; kwargs...)
# LatticeHopsSparse(args...; kwargs...) = LatticeHopsSparse{SparseMatrixCSC{ComplexF64}}(args...; kwargs...)
# LatticeHopsDense{T}(args...; kwargs...) where T<:Matrix{<:Complex} = Dict{Vector{Int}, T}(args...; kwargs...)
# LatticeHopsDense(args...; kwargs...) = LatticeHopsDense{Matrix{ComplexF64}}(args...; kwargs...)

hopdim(hops::LatticeHops) = size(first(values(hops)),1)

empty_hops() = Dict{Vector{Int},SparseMatrixCSC{ComplexF64}}()

function Base.kron(a, b::LatticeHops)
    for (δL, t) in b
        b[δL] = kron(a, t) # add spin degree of freedom # we made the choice to group the matrix in spin blocks
    end

    b
end

function Base.kron(a::LatticeHops, b)
    for (δL, t) in a
        a[δL] = kron(t, b) # add spin degree of freedom # we made the choice to group the matrix in spin blocks
    end

    a
end

function add_hoppings!(hops::LatticeHops, newhops::LatticeHops)
    for (δR, hop)=newhops
        if !haskey(hops,δR)
            hops[δR] = hop
        else
            hops[δR] += hop
        end
    end
    nothing
end

function add_hoppings(hops::LatticeHops, newhops::LatticeHops)
    hops2 = deepcopy(hops)

    add_hoppings!(hops2, newhops)
    hops2
end

function extend_space(hoppings, mode=:nospin) #::LatticeHops
    if mode==:nospin || mode==:id
        return hoppings
    elseif mode==:spinhalf || mode==:σ0
        hoppings = kron(hoppings, σ0)
    elseif mode==:σx
        hoppings = kron(hoppings, σX)
    else
        error("Do not recognize mode '$mode' in extend_space(...).")
    end

    hoppings
end

get_dense(hops) = Dict(δL => Matrix(t) for (δL, t) in hops)

function decide_type(hops::LatticeHops, format)

    if format==:auto
        N = size(first(values(hops)), 1)
        if N < 301
            format=:dense
        end
    end

    if format==:dense
        hops = get_dense(hops)
    end

    hops
end

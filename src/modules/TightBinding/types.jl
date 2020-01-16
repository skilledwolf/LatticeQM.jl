
const Hop  = Pair{Vector{Int}, <:AbstractMatrix}
const Hops = Dict{Vector{Int}, AbstractMatrix}
const AnyHops = Dict{Vector{Int}, <:AbstractMatrix}


DenseHops(kv::Hop...) = Hops(k=>Matrix(v) for (k,v) in kv)
DenseHops(d::AnyHops) = DenseHops(d...)

SparseHops(kv::Hop...) = Hops(k=>sparse(v) for (k,v) in kv)
SparseHops(d::AnyHops) = SparseHops(d...)

hopdim(hops::AnyHops) = size(first(values(hops)),1)

addhops!(hops::AnyHops, newhops::AnyHops...) = merge!(+, hops, newhops...)
addhops(hops::AnyHops, newhops::AnyHops...) = merge(+, hops, newhops...)

function Base.kron(a, b::AnyHops)
    for (δL, t) in b
        b[δL] = kron(a, t) # add spin degree of freedom # we made the choice to group the matrix in spin blocks
    end

    b
end

function Base.kron(a::AnyHops, b)
    for (δL, t) in a
        a[δL] = kron(t, b) # add spin degree of freedom # we made the choice to group the matrix in spin blocks
    end

    a
end

function addspin(hoppings, mode=:nospin) #::AbstractHops
    if mode==:nospin || mode==:id
        return hoppings
    elseif mode==:spinhalf || mode==:σ0
        hoppings = kron(hoppings, σ0)
    elseif mode==:σx
        hoppings = kron(hoppings, σX)
    else
        error("Do not recognize mode '$mode' in addspin(...).")
    end

    hoppings
end

const maximum_dense_size = 300

function decidetype(hops::AnyHops, format)

    if format==:auto
        N = size(first(values(hops)), 1)
        if N < maximum_dense_size + 1
            format=:dense
        end
    end

    if format==:dense
        hops = DenseHops(hops)
    elseif format==:sparse
        hops = SparseHops(hops)
    end

    hops
end

###################################################################################################
# Backwards compatibility
###################################################################################################
export extend_space
@legacyalias addspin extend_space

export decide_type
@legacyalias decidetype decide_type

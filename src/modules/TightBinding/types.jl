
const Hop  = Pair{Vector{Int}, <:AbstractMatrix}
const Hops = Dict{Vector{Int}, AbstractMatrix}
const AnyHops = Dict{Vector{Int}, <:AbstractMatrix}


DenseHops(kv::Hop...) = Hops(k=>Matrix(v) for (k,v) in kv)
DenseHops(d::AnyHops) = DenseHops(d...)

SparseHops(kv::Hop...) = Hops(k=>sparse(v) for (k,v) in kv)
SparseHops(d::AnyHops) = SparseHops(d...)

Hops(M::AbstractMatrix,d::Int=2) = Hops(zeros(Int,d)=>M)

zerokey(h::AnyHops) = zero(first(keys(h)))
getzero(h::AnyHops) = h[zerokey(h)]

hopdim(hops::AnyHops) = size(first(values(hops)),1)

Base.:+(h1::AnyHops, h2::AnyHops) = addhops(h1,h2)
addhops!(hops::AnyHops, newhops::AnyHops...) = merge!(+, hops, newhops...)
addhops(hops::AnyHops, newhops::AnyHops...) = merge(+, hops, newhops...)

Base.:*(h1::AnyHops, h2::AnyHops) = multiplyhops(h1,h2)
Base.:*(h1::AnyHops, h2::AbstractMatrix) = multiplyhops(h1,h2)
Base.:*(h1::AbstractMatrix, h2::AnyHops) = multiplyhops(h1,h2)
multiplyhops(h1::AbstractMatrix, h2::AnyHops) = multiplyhops(Hops(h1),h2)
multiplyhops(h1::AnyHops, h2::AbstractMatrix) = multiplyhops(h1,Hops(h2))
multiplyhops(h1::AnyHops, h2::AnyHops) = Hops(k=>h1[k]*h2[k] for k=intersect(keys(h1),keys(h2)))

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


# const Hop  = Pair{Vector{Int}, T} where T<:AbstractMatrix
# const Hops = Dict{Vector{Int}, AbstractMatrix}
# const AnyHops = Dict{Vector{Int}, T} where T<:AbstractMatrix

const HopVector = Vector{Int}

const Hop{T<:AbstractMatrix}  = Pair{HopVector, T}
const AnyHop    = Hop{T} where {T<:AbstractMatrix}
const DenseHop  = Hop{Matrix}
const SparseHop = Hop{SparseMatrixCSC}

const Hops{T<:AbstractMatrix} = Dict{HopVector, T}
const AnyHops = Hops{T} where {T<:AbstractMatrix}
const DenseHops = Hops{Matrix}
const SparseHops = Hops{SparseMatrixCSC}


# abstract type AbstractHops{T<:AbstractMatrix} end
# abstract type AbstractHamiltonian{T} end

# struct Hamiltonian{T} <: AbstractHamiltonian{T}
#     H::T 
#     mu::Float64
# end

function zerolike(h::AnyHops; format=:auto)
    ρ = Hops{AbstractMatrix}()

    if format==:dense 
        for δL=keys(h)
            ρ[δL] = zeros(ComplexF64, size(h[δL]))
        end
    elseif format==:sparse
        for δL=keys(h)
            ρ[δL] = spzeros(ComplexF64, size(h[δL]))
        end
    else
        for δL=keys(h)
            ρ[δL] = zero(h[δL])
        end
    end

    ρ
end

DenseHops(kv::Hop...) = Hops{Matrix{ComplexF64}}(k=>Matrix(v) for (k,v) in kv)
DenseHops(d::AnyHops) = DenseHops(d...)

SparseHops(kv::Hop...) = Hops{SparseMatrixCSC{Complex{Float64},Int64}}(k=>sparse(v) for (k,v) in kv)
SparseHops(d::AnyHops) = SparseHops(d...)

Hops() = Hops{AbstractMatrix}()
Hops(kv::Hop...) = Hops{AbstractMatrix}(k=>v for (k,v) in kv)
Hops(d::AnyHops) = Hops(d...)
Hops(M::AbstractMatrix,d::Int=2) = Hops(zeros(Int,d)=>M)

zerokey(h::AnyHops) = zero(first(keys(h)))
getzero(h::AnyHops) = h[zerokey(h)]
setzero!(h::AnyHops, M::AbstractMatrix) = (h[zerokey(h)].=M; h)

hopdim(hops::AnyHops) = size(first(values(hops)),1)

Base.:+(h1::AnyHops, h2::AnyHops) = addhops(h1,h2)
addhops!(hops::AnyHops, newhops::AnyHops...) = merge!(+, hops, newhops...)
addhops(hops::AnyHops, newhops::AnyHops...) = merge(+, hops, newhops...)

Base.:*(h::AnyHops, s::Number) = multiplyhops(h,s)
Base.:*(s::Number, h::AnyHops) = multiplyhops(h,s)
Base.:*(h1::AnyHops, h2::AnyHops) = multiplyhops(h1,h2)
Base.:*(h1::AnyHops, h2::AbstractMatrix) = multiplyhops(h1,h2)
Base.:*(h1::AbstractMatrix, h2::AnyHops) = multiplyhops(h1,h2)
multiplyhops(h1::AbstractMatrix, h2::AnyHops) = multiplyhops(Hops(h1),h2)
multiplyhops(h1::AnyHops, h2::AbstractMatrix) = multiplyhops(h1,Hops(h2))
multiplyhops(h1::AnyHops, h2::AnyHops) = Hops(k=>h1[k]*h2[k] for k=intersect(keys(h1),keys(h2)))
multiplyhops(h::AnyHops, s::Number) = Hops(k=>h[k]*s for k=keys(h))
multiplyhops(s::Number, h::AnyHops) = Hops(k=>h[k]*s for k=keys(h))

"""
Naive implementation of combining the linear spaces of two hopping models.
"""
function addhopspace(h1::AnyHops, h2::AnyHops)
    d1 = hopdim(h1)
    d2 = hopdim(h2)
    D = d1+d2

    for δL in keys(h1)
        tmp = spzeros(D,D)
        tmp[1:d1,1:d1] .= h1[δL]
        h1[δL] = tmp
    end
    for δL in keys(h2)
        tmp = spzeros(D,D)
        tmp[d1+1:D,d1+1:D] .= h2[δL]
        h2[δL] = tmp
    end

    addhops(h1,h2)
end

function Base.kron(a, b::AnyHops)
    b2 = deepcopy(b)

    for (δL, t) in b
        b2[δL] = kron(a, t) # add spin degree of freedom # we made the choice to group the matrix in spin blocks
    end

    b2
end

function Base.kron(a::AnyHops, b)
    a2 = deepcopy(a)

    for (δL, t) in a
        a2[δL] = kron(t, b) # add spin degree of freedom # we made the choice to group the matrix in spin blocks
    end

    a2
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

const MAX_DENSE = 500
const MAX_DIAGS = 100

decidetype(hops::AnyHops) = decidetype(hopdim(hops))

function decidetype(N::Int)
    if N < MAX_DENSE + 1
        return :dense
    else
        return :sparse
    end
end


function ensuretype(hops::AnyHops, format=:auto)
    if format==:auto
        # format = decidetype(hops) # old behaviour
        return hops
    end 

    if format==:dense
        hops = DenseHops(hops)
    elseif format==:sparse
        hops = SparseHops(hops)
    end

    hops
end





# function efficientformat(ρ::AnyHops)
#     L = length(̢ρ)
#     @assert L > 0 "Must have at least one hopping element."

#     dims = size(first(values(ρ)))
    
#     A = Array{eltype(valtype(ρ))}(undef, dims..., L)
    
#     keylist = []
#     for (i,δL) in enumerate(keys(ρ))
#         A[:,:,i] .= ρ[δL][:,:]
#         append!(keylist, [δL])
#     end
    
#     A, keylist
# end

# function efficientzero(ρ::AnyHops)
#     L = length(̢ρ)
#     @assert L > 0 "Must have at least one hopping element."

#     dims = size(first(values(ρ)))
    
#     A = zeros(eltype(valtype(ρ)), dims..., L)

#     A, collect(keys(ρ))
# end

# function flexibleformat(A::AbstractArray, keylist::AbstractVector)
#     Dict(L=>Matrix(m) for (L,m)=zip(keylist,eachslice(A; dims=3)))
# end

# function flexibleformat!(ρ::AnyHops, A::AbstractArray, keylist::AbstractVector)
#     for (j_,L)=enumerate(keylist)
#         # ρ[L][:,:] .= m[:,:]
#         # copyto!(ρ[L][:,:], A[:,:,j_])
#         ρ[L][:,:] .= A[:,:,j_]
#     end
#     ρ
# end
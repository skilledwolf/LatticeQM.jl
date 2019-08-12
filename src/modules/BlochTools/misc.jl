
# Define proper iterators for each input type
const kIterable = Union{DiscretePath, <:AbstractMatrix{Float64}, <:AbstractVector{T1}} where {T1<:AbstractVector{Float64}}
eachpoint(kPoints::DiscretePath) = eachcol(kPoints.points)
eachpoint(ks::T) where {T<:AbstractMatrix{Float64}} = eachcol(ks)
eachpoint(ks::T2) where {T1<:AbstractVector{Float64},T2<:AbstractVector{T1}} = ks
points(kPoints::DiscretePath) = kPoints.points
points(ks::T) where {T<:AbstractMatrix{Float64}} = ks
points(ks::T2) where {T1<:AbstractVector{Float64},T2<:AbstractVector{T1}} = matrixcollect(ks)


function sumk(f_k::Function, ks::kIterable)
    ks = points(ks)

    sum(f_k(k) for k=eachcol(ks))/size(ks,2)
end

heaviside(x::AbstractFloat) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x,0.5)))

fermidirac_T(ϵ::AbstractFloat; T::AbstractFloat=0.01) = 1.0/(exp(ϵ/T)+1)
fermidirac_0(ϵ::AbstractFloat) = heaviside(-ϵ)

function fermidirac(ϵ::AbstractFloat; T::AbstractFloat=0.01)
    if T==0.0
        x = fermidirac_0(ϵ)
    else
        x = fermidirac_T(ϵ; T=T)
    end

    x
end

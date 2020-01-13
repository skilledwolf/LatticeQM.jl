using ..Algebra: eigmin_sparse, eigmax_sparse

function tounitrange!(H::T; ϵ::Float64=0.01) where T<:AbstractMatrix
"""
    Scale the hermitian matrix H such that its (bounded) spectrum is [-1.0,1.0].
    ϵ is just to avoid issues at the boundaries of the interval.
"""
    @assert ishermitian(H)

    eigmin0 = issparse(H) ? eigmin_sparse : eigmin
    eigmax0 = issparse(H) ? eigmax_sparse : eigmax

    ϵmin = eigmin0(H)
    ϵmax = eigmax0(H)

    a = (ϵmax-ϵmin)/(2.0-ϵ)
    b = (ϵmax+ϵmin)/2.0

    H[:] = ((H .- b)./a)[:]

    a, b
end

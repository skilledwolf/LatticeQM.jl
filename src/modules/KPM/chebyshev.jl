###############################################################################
# Chebyshev polynomials
###############################################################################

function Tn(n::S, x::T) where {T<:AbstractFloat,S<:Integer}
"""
    Chebyshev polynomial of order n evaluated at x.
"""

  polynomial    = Array{T}(undef, (1,n+1))
  polynomial[1] = one(T)

  for i = 2:n+1
    if i == 2
      polynomial[i] = x
    else
      polynomial[i] = 2*x*polynomial[i-1]-polynomial[i-2]
    end
  end

  return polynomial
end

function Tn(n::S, x::Array{T,1}) where {T<:AbstractFloat,S<:Integer}
"""
    Chebyshev polynomial of order n evaluated for values x=(x_1, x_2, ...).
"""

  polynomial      = Array{T}(undef, (length(x), n+1))
  polynomial[:,1] = ones(T,length(x))

  for i = 2:n+1
    for j = 1:length(x)
      if i == 2
        polynomial[j,i] = x[j]
      else
        polynomial[j,i] = 2*x[j]*polynomial[j,i-1]-polynomial[j,i-2]
      end
    end
  end

  return polynomial
end

###############################################################################
###############################################################################

function evaluate_expansion(μ::T1, x::Vector{T2}) where {T1<:AbstractVector, T2<:AbstractFloat}
"""
    Evaluate the Chebyshev expansion for given expansion coefficients
        μ = (μ_0, μ_1, ...)
    at values x = (x_1, x_2, x_3, ...).

    Parameters a, b scale the function back from unit range into the original interval.
"""
#     x1 = (x.-b)./a
    (1.0./(π .* sqrt.(1.0.-x.^2))) .* (μ[1] .+ 2.0 .* Tn(length(μ)-1,x)[:,2:end] * μ[2:end])
end

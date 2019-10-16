###############################################################################
# Integral kernels to provide optimal cutoff in the Chebyshev expansion
# [see section II.C), Weisse et al, Rev. Mod. Phys. 78 275]
###############################################################################

jackson(j,n) = ( (n-j+1) * cos(π*j/(n+1)) + sin(π*j/(n+1)) * cot(π/(n+1)) )/(n+1)
function jackson!(μ::AbstractVector)
"""
    Apply Jackson Kernel to the Chebychev expansion.
"""
    n = length(μ)
    for j=0:n-1
        μ[j+1] *= jackson(j,n)
    end
    nothing
end


lorentz(j,n; λ) = sinh(λ*(1-j/n))/sinh(λ)
function lorentz!(μ::AbstractVector; kwargs...)
"""
    Apply Lorentz Kernel to the Chebychev expansion.
"""
    n = length(μ)
    for j=0:n-1
        μ[j+1] *= lorentz(j,n; kwargs...)
    end
    nothing
end


using ..Algebra: geteigvals, geteigvecs, geteigen

using ..TightBinding: AnyHops, getbloch

ensureH(H::AbstractMatrix) = H
ensureH(H::Function) = H
ensureH(H::AnyHops) = getbloch(H)


energies(H, args...; kwargs...) = geteigvals(ensureH(H), args...; kwargs...)
wavefunctions(H, args...; kwargs...) = geteigvecs(ensureH(H),args...; kwargs...)
spectrum(H, args...; kwargs...) = geteigen(ensureH(H),args...; kwargs...)



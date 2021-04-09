
import ..Algebra

import ..TightBinding: Hops, getbloch

ensureH(H::AbstractMatrix) = H
ensureH(H::Function) = H
ensureH(H::Hops) = getbloch(H)


energies(H, args...; kwargs...) = Algebra.geteigvals(ensureH(H), args...; kwargs...)
wavefunctions(H, args...; kwargs...) = Algebra.geteigvecs(ensureH(H),args...; kwargs...)
spectrum(H, args...; kwargs...) = Algebra.geteigen(ensureH(H),args...; kwargs...)



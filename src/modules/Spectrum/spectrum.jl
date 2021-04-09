
import ..Algebra

energies(H, args...; kwargs...) = Algebra.geteigvals(H, args...; kwargs...)
wavefunctions(H, args...; kwargs...) = Algebra.geteigvecs(H,args...; kwargs...)
spectrum(H, args...; kwargs...) = Algebra.geteigen(H,args...; kwargs...)

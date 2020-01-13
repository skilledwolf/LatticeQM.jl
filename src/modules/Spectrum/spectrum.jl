
using ..Algebra

energies(args...; kwargs...) = eigvals(args...; kwargs...)
wavefunctions(args...; kwargs...) = eigvecs(args...; kwargs...)
spectrum(args...; kwargs...) = eigen(args...; kwargs...)


###################################################################################################
# Backwards compatibility
###################################################################################################
# export ϵs_dense
# function ϵs_dense(args...; kwargs...)
#     @warn("ϵs_dense(...) was renamed to eigvals_dense(...) and is marked for removal.")
#     Algebra.eigvals_dense(args...; kwargs...)
# end
# export ϵs_sparse
# function ϵs_sparse(args...; kwargs...)
#     @warn("ϵs_sparse(...) was renamed to eigvals_sparse(...) and is marked for removal.")
#     Algebra.eigvals_sparse(args...; kwargs...)
# end
# export Us_dense
# function Us_dense(args...; kwargs...)
#     @warn("Us_dense(...) was renamed to eigvecs_dense(...) and is marked for removal.")
#     Algebra.eigvecs_dense(args...; kwargs...)
# end
# export Us_sparse
# function Us_sparse(args...; kwargs...)
#     @warn("Us_sparse(...) was renamed to eigvecs_sparse(...) and is marked for removal.")
#     Algebra.eigvecs_sparse(args...; kwargs...)
# end
# export Us
# function Us(args...; kwargs...)
#     @warn("Us(...) was renamed to eigvecs(...) and is marked for removal.")
#     Algebra.eigvecs(args...; kwargs...)
# end


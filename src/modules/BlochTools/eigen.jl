"""
    Additional low-level interfaces (acting on arrays)
    ==================================================
"""

eigmax_sparse(H::AbstractMatrix) = eigs(H; nev=1, which=:LR)[1][1] |> real
eigmin_sparse(H::AbstractMatrix) = eigs(H; nev=1, which=:SR)[1][1] |> real
function eigen_sparse(M::AbstractMatrix; nev::Int, sigma::Float64=1e-8, which=:LM, kwargs...)
    eigs(M; nev=nev, sigma=sigma, which=which, kwargs...)
end
function eigvals_sparse(M::AbstractMatrix; nev::Int, sigma::Float64=1e-8, which=:LM, kwargs...)
    eigsol = eigen_sparse(M; nev=nev, sigma=sigma, which=which, kwargs...)
    real.(eigsol[1])
end
function eigvecs_sparse(M::AbstractMatrix; nev::Int, sigma::Float64=1e-8, which=:LM, kwargs...)
    eigsol = eigen_sparse(M; nev=nev, sigma=sigma, which=which, kwargs...)
    real.(eigsol[2])
end

"""
    Interfaces for single k-points
    ==============================
"""

# Dense methods
ϵs_dense(h::Function) = k::AbstractVector{Float64} -> eigvals(Matrix(h(k)))
Us_dense(h::Function, k::AbstractVector{Float64}) = k::AbstractVector{Float64} -> eigvecs(Matrix(h(k)))
function eigen_dense(h::Function)
    function atK(k::AbstractVector{Float64})
        F = eigen(Matrix(h(k)))
        F.values, F.vectors
    end
    atK
end

# Sparse methods
ϵs_sparse(h::Function; kwargs...) = k::AbstractVector{Float64} -> eigvals_sparse(h(k); kwargs...)
Us_sparse(h::Function; kwargs...) = k::AbstractVector{Float64} -> eigvecs_sparse(h(k); kwargs...)
eigen_sparse(h::Function; kwargs...) = k::AbstractVector{Float64} -> eigen_sparse(h(k); kwargs...)

function spectrum(h::Function; format=:dense, kwargs...)
    if format==:dense
        eigen_dense(h; kwargs...)
    elseif format==:sparse
        eigen_sparse(h; kwargs...)
    end
end

################################################################################
################################################################################
"""
    Interfaces for discretized k-path
    =================================
"""
# Define proper iterators for each input type
const kIterable = Union{DiscretePath, <:AbstractMatrix{Float64}, <:AbstractVector{T1}} where {T1<:AbstractVector{Float64}}
eachpoint(kPoints::DiscretePath) = eachcol(kPoints.points)
eachpoint(ks::T) where {T<:AbstractMatrix{Float64}} = eachcol(ks)
eachpoint(ks::T2) where {T1<:AbstractVector{Float64},T2<:AbstractVector{T1}} = ks

energies(h::Function, ks::kIterable)      = Base.Generator(ϵs_dense(h), eachpoint(ks))
wavefunctions(h::Function, ks::kIterable) = Base.Generator(Us_dense(h), eachpoint(ks))
energies_wfs(h::Function, ks::kIterable)  = Base.Generator(eigen_dense(h), eachpoint(ks))

################################################################################
################################################################################
"""
    Calculate and store bands for kpoints ks
    ========================================
"""

using RecursiveArrayTools

matrixcollect(it) = convert(Array, VectorOfArray(collect(it)))
bandmatrix(h::Function, ks::kIterable) = matrixcollect(energies(h, ks))

using Distributed
function bandmatrix_parallel(h::Function, ks)
    matrixcollect(pmap(ϵs_dense(h), eachpoint(ks)))
end

###################################################################################################
###################################################################################################

# function LinearAlgebra.:eigvals(h::Function, point::T; mode=:dense, kwargs...) where {T<:AbstractVector}
#
#     if mode==:dense
#         f = eigvals
#     elseif mode==:sparse
#         f = eigvals_sparse
#     else
#         error("Unknown mode.")
#     end
#
# #     bands = @showprogress 0.1 "Diagonalization" map(x->f(Matrix(h(Vector(x))), kwargs...), eachpoint(points))
# #
# #     Matrix{Float64}(hcat(bands...))
#     f(Matrix(h(Vector(point))), kwargs...)
# end

###################################################################################################
###################################################################################################

function chemical_potential(hamiltonian::Function, ks::AbstractMatrix{Float64}, filling::Float64; type=:sparse)
    if type==:sparse
        out = chemical_potential_sparse_parallel(hamiltonian, ks, filling)
    else
        out = chemical_potential_dense(hamiltonian, ks, filling)
    end
    out
end

function chemical_potential_sparse_parallel(hamiltonian::Function, ks::AbstractMatrix{Float64}, filling::Float64)
"""
    Possible issues:
    - bands can be huge for twisted systems, making the storage of "bands" inefficient.
"""

    min = minimum( pmap(k->eigmin_sparse(hamiltonian(Vector(k))), eachcol(ks)) )
    max = maximum( pmap(k->eigmax_sparse(hamiltonian(Vector(k))), eachcol(ks)) )

    min + filling * (max-min)
end

function chemical_potential_sparse(hamiltonian::Function, ks::AbstractMatrix{Float64}, filling::Float64)
"""
    Possible issues:
    - bands can be huge for twisted systems, making the storage of "bands" inefficient.
"""
    min = minimum(eigmin_sparse(hamiltonian(Vector(k))) for k in eachcol(ks))
    max = maximum(eigmax_sparse(hamiltonian(Vector(k))) for k in eachcol(ks))

    min + filling * (max-min)
end

function chemical_potential_dense(hamiltonian::Function, ks::AbstractMatrix{Float64}, filling::Float64)
"""
    Possible issues:
    - bands can be huge for twisted systems, making the storage of "bands" inefficient.
    - also note that this method performs a full diagonalization on all k points,
      which can become very expensive for large systems
"""

    min = minimum(eigmin(hamiltonian(Vector(k))) for k in eachcol(ks))
    max = maximum(eigmax(hamiltonian(Vector(k))) for k in eachcol(ks))

    min + filling * (max-min)
end

###################################################################################################
###################################################################################################

function bandgap_filling_dense(hamiltonian::Function, ks::AbstractMatrix{Float64}, filling::Float64)

    # Calculate the gap around in which the Fermi level lies
    bands = bandmatrix(hamiltonian, ks)

    ϵmin = minimum(bands)
    ϵmax = maximum(bands)
    ϵfermi = ϵmin + filling * (ϵmax-ϵmin)

    ϵlower = maximum(bands[bands.<= ϵfermi])
    ϵupper = minimum(bands[bands.>= ϵfermi])
    ϵgap = ϵupper - ϵlower

    ϵgap
end

function bandgap_μ_dense(hamiltonian::Function, ks::AbstractMatrix{Float64}, μ::Float64)

    # Calculate the gap around in which the Fermi level lies
    bands = bandmatrix(hamiltonian, ks)

    ϵlower = maximum(bands[bands.<= μ])
    ϵupper = minimum(bands[bands.>= μ])
    ϵgap = ϵupper - ϵlower

    ϵgap
end

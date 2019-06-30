__precompile__()

module LatticeTools
################################################################################
################################################################################

using Base, LinearAlgebra
using RecursiveArrayTools

################################################################################
################################################################################

mutable struct Lattice{T<:AbstractMatrix{Float64}}

    """
        D::Int :                                   Dimension of the lattice
        A::Matrix(D, D) :                          Save primitive lattice vectors as columns
        atoms::Matrix(D, N) :                      atom positions in units of cols(A) "fractional coordinates"

        d::Int :                                   Perpendicular spatial dimensions (total spatial dimensions: d+D)
        atoms_auxspace::Matrix{Float}(M,N):        position of each atom in auxiliary space of dimension M
    """

    A::T # d × d
    atoms::T # d × N # fractional coordinates w.r.t. A

    atoms_aux::T # D × N
    extradimensions::Dict{String,Int}
end

Lattice(A::T) where {T<:AbstractMatrix{Float64}} = Lattice{T}(A, zeros(Float64, size(A)[1], 1))

function Lattice(A::T, atoms::T; extradimensions::Vector{String}=Vector{String}()) where {T<:AbstractMatrix{Float64}}
    @assert size(atoms)[1] == size(A)[1]
    @assert size(A)[1] == size(A)[2]

    L = size(extradimensions)

    atoms_aux = zeros(Float64, L, size(atoms)[2]) # by default there is no perpendicular space
    extradimensions_dict = Dict{String,Int}(key=>index for (index,key) in enumerate(extradimensions))

    Lattice{T}(A, atoms, atoms_aux, extradimensions_dict)
end

function Lattice(A::T, atoms::T, atoms_aux::T; extradimensions::Vector{String}=Vector{String}()) where {T<:AbstractMatrix{Float64}}
    @assert size(atoms)[1] == size(A)[1]
    @assert size(A)[1] == size(A)[2]
    @assert size(atoms_aux)[2] == size(atoms)[2]

    extradimensions_dict = Dict{String,UInt}(key=>index for (index,key) in enumerate(extradimensions))

    Lattice{T}(A, atoms, atoms_aux, extradimensions_dict)
end

################################################################################
################################################################################

include("LatticeTools/lattice_methods.jl")
include("LatticeTools/supercell_methods.jl")

################################################################################
################################################################################
end

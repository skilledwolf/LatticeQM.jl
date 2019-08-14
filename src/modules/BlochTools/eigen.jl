
    # Additional low-level interfaces (acting on arrays)
    # ==================================================

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

    # Interfaces for single k-points
    # ==============================

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

function ϵs(h::Function; format=:dense, kwargs...)
    if format==:dense
        ϵs_dense(h; kwargs...)
    elseif format==:sparse
        ϵs_sparse(h; kwargs...)
    end
end
function spectrum(h::Function; format=:dense, kwargs...)
    if format==:dense
        eigen_dense(h; kwargs...)
    elseif format==:sparse
        eigen_sparse(h; kwargs...)
    end
end

################################################################################
################################################################################

    # Interfaces for discretized k-path
    # ==============================


using RecursiveArrayTools
matrixcollect(it) = convert(Array, VectorOfArray(collect(it)))

energies(h::Function, ks::kIterable)      = Base.Generator(ϵs_dense(h), eachpoint(ks))
wavefunctions(h::Function, ks::kIterable) = Base.Generator(Us_dense(h), eachpoint(ks))
energies_wfs(h::Function, ks::kIterable)  = Base.Generator(eigen_dense(h), eachpoint(ks))


################################################################################
################################################################################

    # Calculate and store bands for kpoints ks
    # ========================================

using Distributed

bandmatrix(h::Function, ks::kIterable) = matrixcollect(energies(h, ks))

function get_bands(h::Function, ks::kIterable; projector=nothing, kwargs...)
    """
        h(k): returns hermitian Hamiltonian at k-point
        ks: collection of k points (see kIterable)
        projector: function that returns a real value for a given (k,ψ,E_k).
            can also be as Vector of such functions.
    """
    ks = points(ks)
    N = size(ks)[2]         # no. of k points
    M = size(h(ks[:,1]))[1] # no. of bands

    bands = zeros(Float64, M, N)

    if projector != nothing
        if isa(projector, Vector)
            obs = zeros(Float64, M, N, length(projector))
        else
            obs = zeros(Float64, M, N, 1)
            projector = [projector]
        end
    end

    Σ = spectrum(h; kwargs...)
    projector = projector

    # Parallized loop
    bands = convert(SharedArray, bands)
    if projector != nothing
        obs = convert(SharedArray, obs)
    end 

    @sync @distributed for j_=1:N
        k = ks[:,j_]
        ϵs, U = Σ(k)

        bands[:,j_] .= ϵs

        if projector != nothing
            for (i_,ψ)=enumerate(eachcol(U))
                for (n_, proj)=enumerate(projector)
                    obs[i_,j_,n_] = proj(k,ψ,ϵs[i_])
                end
            end
        end
    end

    bands = convert(Array, bands)

    if projector == nothing
        obs = nothing
    else
        obs = convert(Array, obs)
    end

    bands, obs
end

function bandmatrix_parallel(h::Function, ks)
    matrixcollect(pmap(ϵs_dense(h), eachpoint(ks)))
end

function groundstate_sumk(ϵs_k::AbstractVector{Float64}, μ::Float64=0.0; kwargs...)
    tmp = 0.0
    for ϵ in ϵs_k
        if ϵ <= μ
            tmp += ϵ
        end
    end

    tmp
end

function groundstate_energy(ϵs::Function, ks::AbstractMatrix{Float64}, μ::Float64=0.0; kwargs...)
    # Σ = ϵs(hamiltonian; format=format)
    L = size(ks)[2]

    ϵGS = @distributed (+) for j=1:L # @todo: this should be paralellized
        tmp = 0.0
        for ϵ in ϵs(ks[:,j])
            if ϵ <= μ
                tmp += ϵ
            end
        end
        tmp
    end

    ϵGS / L
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

function chemical_potential(hamiltonian::Function, ks::AbstractMatrix{Float64}, filling::Float64; T::AbstractFloat=0.0)
    if T==0.0
        μ = chemical_potential_0(hamiltonian, ks, filling)
    else
        μ = chemical_potential_T(hamiltonian, ks, filling; T=T)
    end

    μ
end

function chemical_potential_0(hamiltonian::Function, ks::AbstractMatrix{Float64}, filling::Float64)

    energies = sort(reshape(bandmatrix_parallel(hamiltonian, ks), :))

    i = floor(Int, filling * size(energies,1)) # fill a fraction of states according to ,,filling''

    (energies[i] + energies[i+1])/2.0
end

using NLsolve

function chemical_potential_T(hamiltonian::Function, ks::AbstractMatrix{Float64}, filling::Float64; T::Float64=0.01)

    d = size(hamiltonian(first(eachcol(ks))), 1)
    N = size(ks,2)

    bands = bandmatrix_parallel(hamiltonian, ks)
    energies = sort(reshape(bands, :))

    N_occ = floor(Int, filling * length(bands)) # fill a fraction of states according to ,,filling''

    # Initial guess for T=0
    μ0 = (energies[N_occ] + energies[N_occ+1])/2.0

    ##########################################
    ###  Solve  n(μ_star) = N_occ for μ_star
    ##########################################

    n(μ::AbstractFloat) = sum(fermidirac(ϵ-μ; T=T) for ϵ=bands)/N

    function δn!(δn::T, μ::T) where {T2<:AbstractFloat, T<:AbstractArray{T2}}
        δn[1] = n(μ[1]) - d * filling
        nothing
    end

    sol = nlsolve(δn!, [μ0])
    @assert converged(sol)
    μ = sol.zero[1]

    μ
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

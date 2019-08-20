
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

function ϵs(h::Function; format=:dense, num_bands=nothing, kwargs...)
    if format==:dense
        ϵs_dense(h; kwargs...)
    elseif format==:sparse

        if num_bands==nothing
            ϵs_sparse(h; kwargs...)
        else
            ϵs_sparse(h; nev=num_bands, kwargs...)
        end
    end
end

function spectrum(h::Function; format=:dense, num_bands=nothing, kwargs...)
    if format==:dense
        eigen_dense(h; kwargs...)
    elseif format==:sparse
        if num_bands==nothing
            eigen_sparse(h; kwargs...)
        else
            eigen_sparse(h; nev=num_bands, kwargs...)
        end
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
using ProgressMeter

# bandmatrix(h::Function, ks::kIterable) = matrixcollect(energies(h, ks))

function get_bands(h::Function, ks::kIterable; projector=nothing, num_bands=nothing, kwargs...)
    """
        h(k): returns hermitian Hamiltonian at k-point
        ks: collection of k points (see kIterable)
        projector: function that returns a real value for a given (k,ψ,E_k).
            can also be as Vector of such functions.
    """
    ks = points(ks)
    N = size(ks)[2]         # no. of k points

    if num_bands != nothing
        M = num_bands
    else
        M = size(h(ks[:,1]))[1] # no. of bands
    end

    bands = zeros(Float64, M, N)

    if projector != nothing
        if isa(projector, Vector)
            obs = zeros(Float64, M, N, length(projector))
        else
            obs = zeros(Float64, M, N, 1)
            projector = [projector]
        end
    end

    Σ = spectrum(h; num_bands=num_bands, kwargs...)
    projector = projector

    # Parallized loop
    bands = convert(SharedArray, bands)
    if projector != nothing
        obs = convert(SharedArray, obs)
    end

    ## EXPERIMENTAL USE OF THE PROGRESS BAR. Spoiler: does not work -.-
    p = Progress(N, 0.1, "Computing bands...")
    channel = RemoteChannel(()->Channel{Bool}(10), 1)

    @sync begin
        # this task prints the progress bar
        @async while take!(channel)
            next!(p)
        end

        # this task does the computation
        @async begin
            @distributed for j_=1:N
                k = ks[:,j_]
                ϵs, U = Σ(k)

                bands[:,j_] .= real.(ϵs)

                if projector != nothing
                    for (i_,ψ)=enumerate(eachcol(U))
                        for (n_, proj)=enumerate(projector)
                            obs[i_,j_,n_] = proj(k,ψ,ϵs[i_])
                        end
                    end
                end
            end
            put!(channel, false) # this tells the printing task to finish
        end
    end
    ## END OF EXPERIMENTAL USE OF THE PROGRESSBAR

    bands = convert(Array, bands)

    if projector == nothing
        obs = nothing
    else
        obs = convert(Array, obs)
    end

    bands, obs
end

# function bandmatrix_parallel(h::Function, ks; kwargs...)
#     matrixcollect(pmap(ϵs(h; kwargs...), eachpoint(ks)))
# end

function bandmatrix(h::Function, ks; num_bands=nothing, kwargs...)
    ks = points(ks)
    L = size(ks)[2]         # no. of k points

    if num_bands != nothing
        M = num_bands
    else
        M = size(h(ks[:,1]))[1] # no. of bands
    end


    bands = convert(SharedArray, zeros(Float64, M, L))

    energies = ϵs(h; num_bands=num_bands, kwargs...)
    @sync @distributed for j_=1:L

        bands[:,j_] .= real.(energies(ks[:,j_]))
    end

    bands
end

@fastmath function groundstate_sumk(ϵs_k::AbstractVector{Float64}, μ::Float64=0.0; kwargs...)
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

# TODO: implement chemical_potential!(bands::AbstractMatrix, ...)
# the repeated allocation of bands is problematic for large systems!

function chemical_potential(hamiltonian::Function, ks::AbstractMatrix{Float64}, filling::Float64; T::AbstractFloat=0.0)

    energies = bandmatrix(hamiltonian, ks)[:]
    nk = size(ks,2)

    chemical_potential!(energies, nk, filling; T=T)
end

function chemical_potential!(energies::AbstractVector{Float64}, nk::Int, filling::Float64; T::AbstractFloat=0.0)

    if T==0.0
        μ = chemical_potential_0!(energies, filling)
    else
        μ = chemical_potential_T!(energies, nk, filling; T=T)
    end

    μ
end

function chemical_potential_0!(energies::AbstractVector{Float64}, filling::Float64)

    i = floor(Int, filling * length(energies)) # fill a fraction of states according to ,,filling''

    e1, e2 = partialsort!(energies, i:i+1)

    # sort!(energies)
    # e1, e2 = energies[i:i+1]

    (e1+e2)/2
end

using NLsolve

function chemical_potential_T!(energies::AbstractVector{Float64}, nk::Int, filling::Float64; T::Float64=0.01)

    d = div(length(energies),nk)

    # Initial guess for T=0
    μ0 = chemical_potential_0!(energies, filling)

    ##########################################
    ###  Solve  n(μ_star) = N_occ for μ_star
    ##########################################

    n(μ::AbstractFloat) = sum(fermidirac(ϵ-μ; T=T) for ϵ=energies)/nk

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

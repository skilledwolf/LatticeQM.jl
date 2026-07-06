import LinearAlgebra: det, norm, I, svd
import SharedArrays: SharedArray, sdata

import LatticeQM.Spectrum
import LatticeQM.Eigen
import LatticeQM.Parallel

@inline _eigvecs_at(H, k; kwargs...) = Eigen.geteigvecs(H(k); kwargs...)

# Common eigenvector-grid loop shared between `statesgrid` and `statesgrid1D`.
# Builds H(k) into a per-task scratch via the Spectrum/bands.jl
# `bloch_buffer` / `bloch!` primitives, runs the eigensolver in-place, and
# slices the requested `bandindices` columns into `out` via the caller's
# `placement(out, j) -> view`. Reuses the same `kspace_foreach!` executor
# pipeline as `bandmatrix`, so callers gain :threaded and :auto support and
# scratch reuse for free.
function _eigvecs_grid!(out, H, ks::AbstractMatrix, bandindices, placement::P;
                        multimode=:auto,
                        executor::Union{Nothing,Parallel.Executor}=nothing,
                        format::Symbol=:dense) where {P}
    exec = executor === nothing ? Parallel.to_executor(multimode) : executor
    if exec isa Parallel.DistributedExec
        H = sanatize_distributed_hamiltonian(H)
    end
    Parallel.configure_blas!(exec; verbose=false)

    Parallel.kspace_foreach!(ks, exec;
        scratch_factory = () -> (Hcache = bloch_buffer(H, ks; format=format),)
    ) do scratch, j, k
        Hk = bloch!(scratch.Hcache, H, k)
        _, U = Eigen.geteigen!(Hk)
        @views placement(out, j) .= U[:, bandindices]
    end
    out
end

"""
    statesgrid(H, NX::Int, NY::Int=0, bandindices::AbstractArray=[]; kwargs...)

Evaluates the eigenvectors on a discretized grid (2D Hamiltonian only!) and
stores the result (preserving the grid information). Useful when plaquette
phases need to be calculated.

`statesgrid[i,j,:,k]` is the `k`-th eigenvector at gridpoint `i,j`.

Accepts `multimode`/`executor` and `format` like `bandmatrix`; defaults to
`:auto` and `:dense`.
"""
function statesgrid(H, NX::Int, NY::Int=0, bandindices::AbstractArray=[];
                    multimode=:auto,
                    executor::Union{Nothing,Parallel.Executor}=nothing,
                    format::Symbol=:dense)
    NY = (NY<1) ? NX : NY

    # 2D k-grid (NX × NY matrix of 2-vectors). Flattened column-major into a
    # 2 × (NX*NY) matrix for `kspace_foreach!`; linear index `j` maps back to
    # `(mod1(j, NX), cld(j, NX))`.
    kgrid = [[x;y] for x=range(0; stop=1, length=NX), y=range(0; stop=1, length=NY)]
    midkgrid = [[1/(2*NX)+x;1/(2*NY)+y] for x=range(0; stop=1-1.0/NX, length=NX-1), y=range(0; stop=1-1.0/NY, length=NY-1)]

    ks = reduce(hcat, vec(kgrid)) .+ 1.34e-8

    M1 = dim(H, ks)
    bandindices = isempty(bandindices) ? collect(1:M1) : bandindices
    M2 = length(bandindices)

    exec = executor === nothing ? Parallel.to_executor(multimode) : executor
    out = exec isa Parallel.DistributedExec ?
              SharedArray{ComplexF64}(NX, NY, M1, M2) :
              zeros(ComplexF64, NX, NY, M1, M2)

    placement = (out, j) -> @views out[mod1(j, NX), cld(j, NX), :, :]
    _eigvecs_grid!(out, H, ks, bandindices, placement;
                   multimode=multimode, executor=exec, format=format)

    statesgrid0 = exec isa Parallel.DistributedExec ? sdata(out) : out
    kgrid, midkgrid, statesgrid0
end

function statesgrid1D(H, NX::Int, bandindices::AbstractArray=[];
                      multimode=:auto,
                      executor::Union{Nothing,Parallel.Executor}=nothing,
                      format::Symbol=:dense)
    kgrid = LinRange(0, 1, NX)
    ks = reshape(collect(Float64, kgrid), 1, NX)

    M1 = dim(H, ks)
    bandindices = isempty(bandindices) ? collect(1:M1) : bandindices
    M2 = length(bandindices)

    exec = executor === nothing ? Parallel.to_executor(multimode) : executor
    out = exec isa Parallel.DistributedExec ?
              SharedArray{ComplexF64}(NX, M1, M2) :
              zeros(ComplexF64, NX, M1, M2)

    placement = (out, j) -> @views out[j, :, :]
    _eigvecs_grid!(out, H, ks, bandindices, placement;
                   multimode=multimode, executor=exec, format=format)

    statesgrid0 = exec isa Parallel.DistributedExec ? sdata(out) : out
    kgrid, statesgrid0
end


function L(State1::T, State2::T) where {T<:AbstractArray{<:Complex,N}} where N
    res = det(State1' * State2)
    return res/abs(res)
end

"""
    plaquettephase(S00, S10, S01, S11)

Calculates the phase the (non-abelian) phase winding around a plaquette.

This method is not meant to be called directly, it is used by `berry(statesgrid)`.
"""
function plaquettephase(S00::T, S10::T, S01::T, S11::T) where {T<:AbstractArray{<:Complex, N}} where N
    real(1.0/(2π*1.0im) * log( L(S00, S10) * L(S10, S11) * L(S01, S11)^(-1) * L(S00,S01)^(-1) ))
end

mymod(i::Int,j::Int) = 1+mod(i-1, j)

"""
    berry(statesgrid)

Goes through each plaquette `(i,j),(i+1,j),(i+1,j+1),(i,j+1)` and calculates the (non-abelian) plaquette phase.
`statesgrid` is a four-dimenional array, containing the discretization information and the occupied states.

You can create a statesgrid with `statesgrid(H, nx, ny, bandindices)` or use the wrapper `berry(H, nx, ny, bandindices)`.
"""
function berry(statesgrid::AbstractArray{<:Complex, 4})

    (n,m) = size(statesgrid)[1:2]
    F = zeros(n-1,m-1)

    for i=1:n-1, j=1:m-1
        S00 = statesgrid[i,  j, :, :]
        S10 = statesgrid[i+1,j, :, :]
        S01 = statesgrid[i, j+1, :, :]
        S11 = statesgrid[i+1, j+1, :, :]
        F[i,j] = plaquettephase(S00, S10, S01, S11)
    end

    F
end

"""
    berry(H, NX::Int, NY::Int=0, bandindices::AbstractArray=[])

Convenience method for `berry(statesgrid)`.
"""
function berry(H, NX::Int, NY::Int=0, bandindices::AbstractArray=[1])
    # geteigvecs(k::Vector) -> Matrix{Complex}
    kgrid, midkgrid, statesgrid0 = statesgrid(H, NX, NY, bandindices)
    midkgrid, berry(statesgrid0)
end

import ..Structure: rotation2D

"""
    berryalongpath(H, kpoints)

Calculate the abelian Berry curvature for each band along a path of discrete k points.
It builds little plaquettes along the path between the kpoints[:,i] and kpoints[:,i+1].
"""
function berryalongpath(H, kpoints;
                        multimode=:auto,
                        executor::Union{Nothing,Parallel.Executor}=nothing,
                        format::Symbol=:dense)

    N = size(kpoints, 2)            # number of k-points along the path
    M = dim(H, kpoints)             # Hilbert dimension

    # Padding endpoints so each interior point has neighbours on both sides
    # for plaquette construction.
    k0 = kpoints[:, 1] - (kpoints[:, 2] - kpoints[:, 1])
    kNp1 = kpoints[:, N] + (kpoints[:, N] - kpoints[:, N - 1])
    kpoints = hcat(k0, kpoints[:, :], kNp1)

    # Random in-plane shift so the path doesn't sit on top of a band touching
    # (Dirac point, etc.) where det(U' * U') has a phase singularity.
    kpoints .+= 0.10 * norm(kpoints[:, 2] - kpoints[:, 1]) * rand(2)

    centers = view(kpoints, :, 2:N+1)

    exec = executor === nothing ? Parallel.to_executor(multimode) : executor
    if exec isa Parallel.DistributedExec
        H = sanatize_distributed_hamiltonian(H)
    end
    Parallel.configure_blas!(exec; verbose=false)

    berryc = exec isa Parallel.DistributedExec ?
                 SharedArray{Float64}(M, N) :
                 zeros(Float64, M, N)

    Parallel.kspace_foreach!(centers, exec;
        scratch_factory = () -> (Hcache = bloch_buffer(H, centers; format=format),)
    ) do scratch, j, _kc_view
        # `j` is the centre index 1..N; the corresponding column in `kpoints`
        # (after padding) is `j+1`. We re-read `kpoints` here to keep the
        # neighbour-aware plaquette geometry that the path-derivative needs.
        kc = kpoints[:, j+1]
        δkR = (kpoints[:, j+2] - kc) / 2
        δkL = (kpoints[:, j]   - kc) / 2
        δkU = (rotation2D( π/2) * δkR + rotation2D(-π/2) * δkL) / 2
        δkD = (rotation2D(-π/2) * δkR + rotation2D( π/2) * δkL) / 2

        # Eigenvectors at the four plaquette corners. The shared `Hcache` is
        # overwritten between corners; `Eigen.geteigen!` returns a fresh `U`
        # matrix per call, so the four `Us` are independent of each other and
        # of the Hcache state once collected.
        Us = map((δkR, δkU, δkD, δkL)) do δk
            Hk = bloch!(scratch.Hcache, H, kc .+ δk)
            _, U = Eigen.geteigen!(Hk)
            U
        end

        for i_ in 1:M
            berryc[i_, j] = plaquettephase((U[:, i_] for U in Us)...)
        end
    end

    exec isa Parallel.DistributedExec ? sdata(berryc) : berryc
end


"""
    getberry!(bands, h, ks)

Calculates and appends the data from `berryalongpath(h,ks)` to the data object `bands`.
This is a convenience method that is useful when plotting band diagrams with Berry curvatures colored it.
"""
function getberry!(bands::BandData, h, ks)
    
    obs = Array(berryalongpath(h, ks))
    obs = reshape(obs, (size(obs)...,1))

    if bands.obs == nothing
        bands.obs = obs
    else
        bands.obs = cat(bands.obs, obs, dims=(3))
    end

    nothing
end

@deprecate NestedWilson2D Wilson2D
@deprecate NestedWilsonWannier2D WilsonWannier2D

function Wilson2D(H, NX::Int, NY::Int=0, bandindices=[1])

    _, _, statesgrid0 = statesgrid(H, NX, NY, bandindices)
    Wilson2D(statesgrid0[1:(NX-1),1:(NY-1),:,:])
end

function Wilson2D(statesgrid::AbstractArray{ComplexF64,4})
    NX, NY, M1, M2 = size(statesgrid)

    en1, U1= WilsonSlice1D(statesgrid, 1)
    en2, U2 = WilsonSlice1D(statesgrid, 2)

    en1, U1, en2, U2
end


function WilsonWannier2D(H, NX::Int, NY::Int=0, bandindices=[1])

    _, _, statesgrid0 = statesgrid(H, NX, NY, bandindices)
    WilsonWannier2D(statesgrid0[1:(NX-1),1:(NY-1),:,:])
end

import Statistics

function WilsonWannier2D(statesgrid::AbstractArray{ComplexF64,4})
    NX, NY, M1, M2 = size(statesgrid)

    newstatesgrid1 = similar(statesgrid)
    for i_=1:NX, j_=1:NY
        @views _, U = Wilson1D(statesgrid[i_,:,:,:], j_-1)
        newstatesgrid1[i_,j_,:,:] = statesgrid[i_,j_, :, :] * U
    end

    newstatesgrid2 = similar(statesgrid)
    for i_=1:NX, j_=1:NY
        @views _, U = Wilson1D(statesgrid[:,j_,:,:], i_-1)
        newstatesgrid2[i_,j_,:,:] = statesgrid[i_,j_, :, :] * U
    end

    @views p1 = [Statistics.mean(Wilson1D(newstatesgrid1[:,j_,:,n_]) for j_=1:NY) for n_=1:M2]
    @views p2 = [Statistics.mean(Wilson1D(newstatesgrid2[i_,:,:,n_]) for i_=1:NX) for n_=1:M2]
    
    p1, p2
end


function Wilson1D(H, NX::Int, bandindices::AbstractArray=[])
    _, statesgrid = statesgrid1D(H, NX, bandindices)

    en, U = Wilson1D(statesgrid)

    en, U, statesgrid
end

function Wilson1D(statesgrid::AbstractArray{ComplexF64,3}, j0::Int=0)
    NY, M1, M2 = size(statesgrid)

    F = Matrix((1.0+0im)*I, (M2,M2))
    for j_=1:NY
        SVD = svd(statesgrid[mymod(j0+j_,NY), :, :]' * statesgrid[mymod(j0+j_+1,NY), :, :])#1+mod(j_-1+1,NY-1)
        F = F * (SVD.U * SVD.Vt)
        # F *= statesgrid[mymod(j0+j_,NY), :, :]' * statesgrid[mymod(j0+j_+1,NY), :, :]
    end

    en, U = Eigen.geteigen(F)
    en = angle.(en)./(2π)
    perm = sortperm(en)

    en[perm], U[:,perm]
end

function Wilson1D(statesgrid::AbstractArray{ComplexF64,2}, j0::Int=0)
    NY, M1= size(statesgrid)

    F = 1.0+0im
    for j_=1:NY
        F = F * (statesgrid[mymod(j0+j_,NY),:]' * statesgrid[mymod(j0+j_+1,NY),:])
    end

    angle(F)/(2π)
end

function WilsonSlice1D(statesgrid::AbstractArray{ComplexF64,4})
    NX, NY, M1, M2 = size(statesgrid)

    en1 = Array{Float64}(undef, (NX,M2))
    U1 = Array{ComplexF64}(undef, (NX,M2,M2))
    # W1 = Array{ComplexF64}(undef, (NX,NY,M1,M2))
    for i_=1:NX
        @views en, U = Wilson1D(statesgrid[i_,  :, :, :])

        en1[i_,:] = en
        U1[i_,:,:] = U
    end

    en1, U1
end

function WilsonSlice1D(statesgrid::AbstractArray{ComplexF64,4}, dim::Int)
    myperm(n::Int, N::Int) = (v = collect(1:N); v[n] = 1; v[1]=n; v)

    if dim == 1
        return @views WilsonSlice1D(statesgrid)
    else
        perm = myperm(dim,ndims(statesgrid))
        @views en, U= WilsonSlice1D(permutedims(statesgrid, perm))
        # W = permutedims(W, perm)
        return en, U#, W
    end
end



"""
    getcherns(wavefunctions::Function, NX::Int, NY::Int=0, bands::AbstractArray=[])

Returns the chern numbers of wavefunctions corresponding to all bands in bands,
where NX and NY denote the coarseness of discretization of k-space.
For bands=[] it returns all chern numbers.
"""
function getcherns(H, NX::Int, NY::Int=0, bands::AbstractArray=[])
    _, _, grid = statesgrid(H, NX, NY, bands)
    cherns = [sum(Spectrum.berry(grid[:,:,:,[i]])) for i=1:size(grid,4)]
    cherns
end


"""
    getwindnum(wavefunctions::Function, NX::Int, NY::Int=0, bandnr::Integer=1)

Returns the winding number (according to Rudner et al. 2013) of the band which is on position bandnr in the spectrum.
NX and NY denote the coarsness of the discretization in k-space and wavefunctions is a function returning the wavefunctions of the problem.
(NB: The winding number of a band is the sum of the chern numbers of all the bands below it including its own chern number.
For undriven systems or systems which are driven but "normal" it is equal to the chern number of said band.
It only makes sense when looking at Floquet bands ie quasi-energies.)
"""
function getwindnum(H, NX::Int, NY::Int=0, bandnr::Integer=1)
    return sum(getcherns(H, NX, NY, collect(1:bandnr)))
end


"""
    spectralgaps(H, NX::Int, NY::Int=NX; gaptol=1e-6, multimode=:auto, executor=nothing, format=:dense)

Global spectral gaps of a two-dimensional Bloch Hamiltonian `H`, found on an
`NX × NY` k-grid from **eigenvalues only** (no eigenvectors) — the cheap
companion to [`gapcherns`](@ref).

A gap below band `n` is recorded when the grid-minimum of band `n+1` exceeds the
grid-maximum of band `n` by more than `gaptol` (in the energy units of `H`).

!!! warning "Grid false gaps"
    A band touching that falls *between* grid points is misreported as a
    finite gap of width ~`v·Δk` (`v` the band velocity, `Δk = 2π/NX`). The
    π-flux (Φ = 1/2) honeycomb Hofstadter Hamiltonian is a concrete case: its
    Dirac touchings show up as ~0.2|t| "gaps" on a 15×15 grid. Such false gaps
    shrink with `NX` and, in a Hofstadter context, typically violate the
    Diophantine constraint (no admissible `(s, C)` exists) — both are useful
    tells.

Returns a `Vector` of `NamedTuple`s `(; n, elo, ehi)`, ordered by energy: `n`
bands lie below a gap spanning `[elo, ehi]`.
"""
function spectralgaps(H, NX::Int, NY::Int=NX; gaptol::Real=1e-6,
                      multimode=:auto,
                      executor::Union{Nothing,Parallel.Executor}=nothing,
                      format::Symbol=:dense)
    # Same k-grid construction as `statesgrid`, so band indices here label the
    # same energy-ordered states a subsequent `statesgrid`/Berry pass returns.
    kgrid = [[x; y] for x in range(0; stop=1, length=NX), y in range(0; stop=1, length=NY)]
    ks = reduce(hcat, vec(kgrid)) .+ 1.34e-8
    D = dim(H, ks)
    E = reshape(bandmatrix(H, ks; multimode=multimode, executor=executor,
                           format=format, hidebar=true)[1], D, NX, NY)
    emin = vec(minimum(E; dims=(2, 3)))     # per-band grid minimum
    emax = vec(maximum(E; dims=(2, 3)))     # per-band grid maximum

    gaps = @NamedTuple{n::Int, elo::Float64, ehi::Float64}[]
    for n in 1:D-1
        emin[n+1] - emax[n] > gaptol || continue
        push!(gaps, (; n=n, elo=emax[n], ehi=emin[n+1]))
    end
    gaps
end

"""
    gapcherns(H, NX::Int, NY::Int=NX; gaptol=1e-6, multimode=:auto, executor=nothing, format=:dense)

Chern numbers of the spectral gaps of a two-dimensional Bloch Hamiltonian `H`.

The gaps are located with [`spectralgaps`](@ref); for each gap the returned
Chern number is the cumulative, **non-abelian** Berry flux of the occupied
manifold (the lowest `n` bands),

    C(n) = sum(berry(states[:, :, :, 1:n])),

i.e. the Hall conductance σxy (in units of e²/h) when the chemical potential
lies in the gap. This manifold quantity is gauge-robust even where individual
bands touch — unlike the single-band Chern numbers from [`getcherns`](@ref),
which are ill-defined for the touching subbands of e.g. a Hofstadter spectrum.

Returns a `Vector` of `NamedTuple`s `(; n, elo, ehi, chern)`, ordered by
energy: `n` bands lie below a gap spanning `[elo, ehi]` with Chern number
`chern`.

The grid must resolve the Berry curvature; for the densely packed magnetic
subbands of a Hofstadter problem `NX ≳ 18` is a safe default (coarser grids can
misreport wide-Chern narrow gaps).

See also [`spectralgaps`](@ref), [`getcherns`](@ref), [`berry`](@ref),
`Operators.hofstadter_cherns`.
"""
function gapcherns(H, NX::Int, NY::Int=NX; gaptol::Real=1e-6,
                   multimode=:auto,
                   executor::Union{Nothing,Parallel.Executor}=nothing,
                   format::Symbol=:dense)
    gaps = spectralgaps(H, NX, NY; gaptol=gaptol, multimode=multimode,
                        executor=executor, format=format)
    out = @NamedTuple{n::Int, elo::Float64, ehi::Float64, chern::Int}[]
    isempty(gaps) && return out

    # Eigenvectors (all bands) for the cumulative Berry flux of each manifold.
    _, _, states = statesgrid(H, NX, NY;
                              multimode=multimode, executor=executor, format=format)
    for g in gaps
        C = round(Int, sum(berry(@view states[:, :, :, 1:g.n])))
        push!(out, (; g.n, g.elo, g.ehi, chern=C))
    end
    out
end


##########################################################################
# Some crude (but working) plotting routines:

# using Statistics
#
# function BZheatmap2(lat, cdata::T) where T<:AbstractMatrix{Float64}
#     NX, NY = size(cdata) .+ 1
#
#     p = plot()
#
#     K(k::T) where T<:AbstractArray{Float64,N} where N = Structure.getB(lat) * k
#
#     points = hcat([ K([(i+0.5)/NX, (j+0.5)/NY]) for i=0:NX-2 for j=0:NY-2]...)
#
#     max = Statistics.quantile(abs.(cdata)[:], 0.95)
#     clim = (-max,max)
#
#     scatter(points[1,:], points[2,:], markersize=3, colorbar=true, markercolor=:RdBu, markerstrokewidth=0, clim=clim, zcolor=cdata[:], legend=:none, seriestype=:heatmap, aspectratio=:equal)
# end
#
# function BZheatmap(lat, cdata::T) where T<:AbstractMatrix{Float64}
#     NX, NY = size(cdata) .+ 1
#
#     p = plot()
#
#     K(k::T) where T<:AbstractArray{Float64,N} where N = Structure.getB(lat) * k
#
#     shapes = Vector{Plots.Shape}()
#
#     for i=0:NX-2, j=0:NY-2
#
#         kpolygon = kparam([ [i/NX; j/NY] [(i+1)/NX; j/NY] [(i+1)/NX; (j+1)/NY] [i/NX; (j+1)/NY]])
#
#         append!(shapes, [Plots.Shape(kpolygon[1,:], kpolygon[2,:])])
#
# #         plot!(p, Plots.Shape(kpolygon[1,:], kpolygon[2,:]), zcolor=cdata[i+1,j+1], fill=(true,cgrad(:grays,[0.,0.1,1.0])))
#     end
#
#     plot!(p, shapes, zcolor=cdata, fill=:RdYlBu, legend=:none)
# end

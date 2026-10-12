import LinearAlgebra: det, norm, I, svd, diag
import SharedArrays: SharedArray, sdata

import LatticeQM.Spectrum
import LatticeQM.Eigen
import LatticeQM.Parallel

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
        _, U = Eigen.geteigen!(Eigen.checked_hermitian(Hk))
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
    # Plaquette centers: grid spacing is 1/(N−1) (`range(0,1,length=N)`), so
    # the (N−1) centers sit at (i+½)/(N−1). The old offset 1/(2N) with
    # spacing (1−1/N)/(N−2) drifted off the true centers (for NX=3 it gave
    # {1/6, 5/6} instead of {1/4, 3/4}), skewing Berry-curvature heatmap axes.
    midkgrid = [[x;y] for x=range(0.5/(NX-1); stop=1-0.5/(NX-1), length=NX-1), y=range(0.5/(NY-1); stop=1-0.5/(NY-1), length=NY-1)]

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

    # Tiny fixed in-plane shift so the path doesn't sit exactly on a band
    # touching (Dirac point at a high-symmetry k, etc.) where det(U'*U) has a
    # phase singularity. Deterministic (results must be reproducible) and
    # ~1e-6 of a path step (the old 10%-of-a-step *random* shift evaluated
    # the curvature visibly off the requested path — systematically wrong
    # wherever the curvature is peaked, which is exactly where it matters).
    kpoints .+= 1.34e-6 * norm(kpoints[:, 2] - kpoints[:, 1]) .* [1.0, 0.7071067811865476]

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
            _, U = Eigen.geteigen!(Eigen.checked_hermitian(Hk))
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

function Wilson2D(H, NX::Int, NY::Int=NX, bandindices=[1]; kwargs...)
    NY = NY < 1 ? NX : NY
    NX > 1 && NY > 1 || throw(ArgumentError("Wilson grids need at least two points per direction"))
    _, _, statesgrid0 = statesgrid(H, NX, NY, bandindices; kwargs...)
    Wilson2D(@view statesgrid0[1:(NX-1),1:(NY-1),:,:])
end

function Wilson2D(statesgrid::AbstractArray{ComplexF64,4})
    NX, NY, M1, M2 = size(statesgrid)

    en1, U1= WilsonSlice1D(statesgrid, 1)
    en2, U2 = WilsonSlice1D(statesgrid, 2)

    en1, U1, en2, U2
end


function WilsonWannier2D(H, NX::Int, NY::Int=NX, bandindices=[1];
                         wannier_cuts=(nothing, nothing), kwargs...)
    NY = NY < 1 ? NX : NY
    NX > 1 && NY > 1 || throw(ArgumentError("Wilson grids need at least two points per direction"))
    _, _, statesgrid0 = statesgrid(H, NX, NY, bandindices; kwargs...)
    WilsonWannier2D(@view(statesgrid0[1:(NX-1),1:(NY-1),:,:]); wannier_cuts)
end

import Statistics

"""
    WilsonWannier2D(statesgrid; wannier_cuts=(nothing, nothing)) -> (p1, p2)

Nested polarizations of gapped Wannier sectors on a periodic, uniformly
sampled grid **without repeated endpoints**. `p1` is along the first grid
direction within sectors of the second-direction Wilson loop; `p2` reverses
the directions. Each entry is a sector's total polarization modulo 1, not
an individual band's polarization inside a degenerate sector.

`wannier_cuts[d]` specifies increasing phase cuts (in units of 2π, spanning
less than one period) bounding the sectors for `p[d]`. For example,
`[-0.25, 0.25]` separates bands near 0 and 1/2, including across ±1/2.
By default, cuts are placed midway between distinct Wilson phases at the
first transverse point; numerically degenerate phases form one sector.
Explicit cuts can group several bands, allowing crossings within a sector.
A single cut selects the entire occupied space.

Sector ranks must stay constant and cuts must avoid the sampled spectrum.
Singular links and unresolved sector gaps throw `ArgumentError`; refine the
grid or choose cuts in persistent Wannier gaps. No finite grid certifies an
unsampled gap. Nested determinant phases are unwrapped before averaging;
a nonzero winding or an ambiguous half-period step precludes this periodic
average and is rejected. The Hamiltonian overload removes repeated endpoints
and accepts the execution keywords of [`statesgrid`](@ref).
"""
function WilsonWannier2D(statesgrid::AbstractArray{ComplexF64,4};
                         wannier_cuts=(nothing, nothing))
    length(wannier_cuts) == 2 || throw(ArgumentError("supply Wannier cuts for both directions"))
    _nestedwilson(statesgrid, wannier_cuts[1]),
        _nestedwilson(PermutedDimsArray(statesgrid, (2, 1, 3, 4)), wannier_cuts[2])
end

function _nestedwilson(statesgrid, cuts)
    NX, NY, _, nb = size(statesgrid)
    NX > 0 && NY > 0 && nb > 0 || throw(ArgumentError("a nested loop needs a nonempty grid and occupied space"))
    frames = similar(statesgrid)
    ranks = Int[]
    for i in 1:NX
        links = _wilsonlinks(@view statesgrid[i, :, :, :])
        centres, U = _wilsonspectrum(links)
        if i == 1
            if cuts === nothing
                c, _ = _wilsonblocks(centres)
                cuts = (c .+ vcat(last(c) - 1, c[1:end-1])) ./ 2
            end
            cuts = collect(Float64, cuts)
            !isempty(cuts) && all(isfinite, cuts) && all(diff(cuts) .> 0) &&
                last(cuts) - first(cuts) < 1 || throw(ArgumentError(
                    "Wannier cuts must increase and span less than one period"))
        end
        if length(cuts) > 1
            minimum(abs(cis(2π*c) - cis(2π*b)) for c in centres, b in cuts) > _WILSON_TOL ||
                throw(ArgumentError("Wannier sector gap is closed or unresolved at a phase cut"))
        end
        labels = [searchsortedlast(cuts .- first(cuts), mod(c - first(cuts), 1)) for c in centres]
        groups = [findall(==(b), labels) for b in eachindex(cuts)]
        if i == 1
            ranks = length.(groups)
            all(>(0), ranks) || throw(ArgumentError("Wannier sectors must be nonempty"))
        else
            length.(groups) == ranks || throw(ArgumentError(
                "Wannier sector dimensions change; choose cuts in persistent gaps or refine the grid"))
        end
        U = U[:, reduce(vcat, groups)]
        for j in 1:NY
            frames[i, j, :, :] = statesgrid[i, j, :, :] * U
            U = links[j]' * U
        end
    end
    edges = cumsum(vcat(0, ranks))
    map(eachindex(ranks)) do b
        bands = edges[b]+1:edges[b+1]
        phases = [sum(first(Wilson1D(@view frames[:, j, :, bands]))) for j in 1:NY]
        _periodic_phase_mean(phases)
    end
end

function _periodic_phase_mean(phases)
    steps = angle.(cis.(2π .* (circshift(phases, -1) .- phases))) ./ (2π)
    any(abs.(abs.(steps) .- 0.5) .<= _WILSON_TOL) && throw(ArgumentError(
        "nested polarization has an unresolved half-period step; refine the grid"))
    abs(sum(steps)) < _WILSON_TOL || throw(ArgumentError(
        "nested polarization winds; a single periodic average is not defined"))
    unwrapped = first(phases) .+ cumsum(vcat(0.0, steps[1:end-1]))
    angle(cis(2π * Statistics.mean(unwrapped))) / (2π)
end


function Wilson1D(H, NX::Int, bandindices::AbstractArray=[])
    _, statesgrid = statesgrid1D(H, NX, bandindices)

    en, U = Wilson1D(statesgrid)

    en, U, statesgrid
end

# Orthonormal frames have overlap singular values in [0, 1]. Below this
# numerical threshold the polar factor cannot define a reliable transport.
const _WILSON_TOL = sqrt(eps(Float64))

function _wilsonlinks(statesgrid::AbstractArray{ComplexF64,3})
    N, _, nb = size(statesgrid)
    N > 0 && nb > 0 || throw(ArgumentError("a Wilson loop needs a nonempty grid and occupied space"))
    map(1:N) do j
        overlap = @views statesgrid[j, :, :]' * statesgrid[mod1(j+1, N), :, :]
        S = svd(overlap)
        minimum(S.S) > _WILSON_TOL || throw(ArgumentError(
            "singular overlap on Wilson link $j; refine the momentum grid or check the selected subspace"))
        S.U * S.Vt
    end
end

function _wilsonspectrum(links, j0::Int=0)
    F = Matrix{ComplexF64}(I, size(first(links), 1), size(first(links), 1))
    for j in eachindex(links)
        F = F * links[mod1(j0+j, length(links))]
    end
    # A unitary matrix is normal. Schur vectors preserve orthonormality even
    # at degeneracy, where generic eigenvectors can be nearly parallel.
    S = schur(F)
    centres = angle.(S.values) ./ (2π)
    p = sortperm(centres)
    centres[p], S.Z[:, p]
end

# Connected clusters on the phase circle, including a cluster crossing ±π.
function _wilsonblocks(centres)
    groups = [Int[1]]
    factors = cis.(2π .* centres)
    for j in 2:length(centres)
        if abs(factors[j] - factors[j-1]) <= _WILSON_TOL
            push!(last(groups), j)
        else
            push!(groups, [j])
        end
    end
    if length(groups) > 1 && abs(first(factors) - last(factors)) <= _WILSON_TOL
        append!(first(groups), pop!(groups))
    end
    c = [angle(sum(factors[g])) / (2π) for g in groups]
    p = sortperm(c)
    c[p], groups[p]
end

"""
    Wilson1D(statesgrid, j0=0)

Wilson loop of orthonormal occupied frames, using SVD-unitarized overlaps.
The matrix form returns `(centres, eigenvectors)` at base point `j0+1`
(cyclically wrapped); the single-band form returns its centre modulo 1.
Numerically singular links (smallest singular value ≤ `sqrt(eps(Float64))`)
throw `ArgumentError`: their transport is undefined or unresolved by the grid.
"""
function Wilson1D(statesgrid::AbstractArray{ComplexF64,3}, j0::Int=0)
    _wilsonspectrum(_wilsonlinks(statesgrid), j0)
end

function Wilson1D(statesgrid::AbstractArray{ComplexF64,2}, j0::Int=0)
    only(first(Wilson1D(reshape(statesgrid, size(statesgrid)..., 1), j0)))
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




################################################################################
# Operator-resolved Wilson loops (holonomy-first and sector-first orderings)
################################################################################

"""
    wilsonweights(statesgrid, O; kgrid=nothing)
        -> (; centres, weights, multiplicities)

Holonomy-first operator weights of a 1D Wilson loop. `statesgrid[j, :, n]`
is the n-th orthonormal occupied state at point j. `O` must be Hermitian in
this fixed, periodic orbital basis.

Each returned entry describes a Wilson eigenspace: its centre modulo 1,
loop-averaged **total** expectation of `O`, and dimension. Wilson eigenvalues
separated by at most `sqrt(eps(Float64))` on the unit circle are grouped,
including across ±π. A rank-one block gives the usual individual band weight;
a degenerate block gives its basis-independent trace, not arbitrary individual
weights. For example, `O = I` gives `weights == multiplicities`.

Diagonalize once and parallel-transport the orthonormal eigenspaces using the
same unitary links as [`Wilson1D`](@ref). Without `kgrid`, samples have equal
weight and should cover one period without a repeated endpoint. Supplying
increasing reduced momenta spanning at most one period uses periodic
trapezoidal quadrature, including correct weights for repeated endpoints.
The Hamiltonian wrapper supplies its grid automatically.

Singular links throw `ArgumentError`. Centres and weights can vary continuously
without protecting symmetry. Weighted centres require a branch convention
when the weights are not integers; no quantized spin/layer invariant is implied.
"""
function wilsonweights(statesgrid::AbstractArray{ComplexF64,3}, O::AbstractMatrix; kgrid=nothing)
    N, norb, nb = size(statesgrid)
    size(O) == (norb, norb) || throw(DimensionMismatch("O must act on the orbital space"))
    O = Eigen.checked_hermitian(O)
    links = _wilsonlinks(statesgrid)
    centres, U = _wilsonspectrum(links)
    centres, groups = _wilsonblocks(centres)
    quadrature = _wilsonquadrature(N, kgrid)
    weights = zeros(nb)
    for j in 1:N
        wb = @views statesgrid[j, :, :] * U
        weights .+= quadrature[j] .* real.(diag(wb' * O * wb))
        U = links[j]' * U
    end
    (; centres, weights=[sum(weights[g]) for g in groups], multiplicities=length.(groups))
end

function _wilsonquadrature(N, kgrid)
    kgrid === nothing && return fill(1/N, N)
    length(kgrid) == N || throw(DimensionMismatch("kgrid must have one momentum per frame"))
    k = collect(Float64, kgrid)
    all(isfinite, k) && all(diff(k) .> 0) && last(k) - first(k) <= 1 ||
        throw(ArgumentError("kgrid must increase and span at most one reduced-momentum period"))
    steps = vcat(diff(k), 1 - (last(k) - first(k)))
    (steps .+ circshift(steps, 1)) ./ 2
end

"""
    sectorwilson(statesgrid::AbstractArray{ComplexF64,3}, O::AbstractMatrix)
        -> ((negative=centres_minus, positive=centres_plus), Ogap)

Sector-first Wilson loops of the negative and positive spectral subspaces of
`V(j)' * O * V(j)` within the occupied space. `O` must be Hermitian and the
occupied frames orthonormal. Each sector has its own non-Abelian Wilson loop,
computed by [`Wilson1D`](@ref); its centres are sorted and defined modulo 1.
An empty sector has an empty vector of centres. Degeneracies and crossings
within a same-sign sector are allowed.

`Ogap` is the minimum absolute projected eigenvalue over the sampled grid.
Throw `ArgumentError` if it is at most `sqrt(eps(Float64)) * norm(O)`, if the
sector dimensions change, or if a sector has a singular overlap link. An open
operator gap defines the sign sectors; it does not by itself quantize their
Wilson phases, nor guarantee that a coarse grid resolves the loop.
"""
function sectorwilson(statesgrid::AbstractArray{ComplexF64,3}, O::AbstractMatrix)
    N, norb, nb = size(statesgrid)
    N > 0 && nb > 0 || throw(ArgumentError("a Wilson loop needs a nonempty grid and occupied space"))
    size(O) == (norb, norb) || throw(DimensionMismatch("O must act on the orbital space"))
    O = Eigen.checked_hermitian(O)
    gap_tol = _WILSON_TOL * norm(O)
    secs = similar(statesgrid)
    Ogap = Inf
    nminus = 0
    for j in 1:N
        V = @view statesgrid[j, :, :]
        ev, U = Eigen.geteigen(Eigen.checked_hermitian(V' * O * V))
        Ogap = min(Ogap, minimum(abs, ev))
        Ogap > gap_tol || throw(ArgumentError("projected operator gap is closed or numerically unresolved"))
        n = count(<(0), ev)
        if j == 1
            nminus = n
        elseif n != nminus
            throw(ArgumentError("operator sector dimensions change along the loop; refine the grid to resolve the gap closure"))
        end
        secs[j, :, :] = V * U
    end
    negative = nminus == 0 ? Float64[] : first(Wilson1D(@view secs[:, :, 1:nminus]))
    positive = nminus == nb ? Float64[] : first(Wilson1D(@view secs[:, :, nminus+1:nb]))
    (negative=negative, positive=positive), Ogap
end

"""
    WilsonWeights1D(H, NX::Int, O::AbstractMatrix, bandindices; kwargs...)

Evaluate [`wilsonweights`](@ref) on explicitly selected eigenstates of the
1D Hamiltonian `H`. Returns `(; centres, weights, multiplicities)`. The grid
includes both endpoints, counted with half weight in the loop average.
`bandindices` must be nonempty; the occupied bands are never inferred.
Accepts the execution keywords of [`statesgrid1D`](@ref).
"""
function WilsonWeights1D(H, NX::Int, O::AbstractMatrix, bandindices::AbstractArray; kwargs...)
    isempty(bandindices) && throw(ArgumentError("select the occupied bandindices explicitly"))
    NX > 1 || throw(ArgumentError("a Wilson grid needs at least two points"))
    kgrid, statesgrid = statesgrid1D(H, NX, bandindices; kwargs...)
    wilsonweights(statesgrid, O; kgrid)
end

"""
    SectorWilson1D(H, NX::Int, O::AbstractMatrix, bandindices; kwargs...)

Evaluate [`sectorwilson`](@ref) on explicitly selected eigenstates of the
1D Hamiltonian `H`. Returns `((negative=centres_minus, positive=centres_plus), Ogap)`.
`bandindices` must be nonempty; the occupied bands are never inferred.
Accepts the execution keywords of [`statesgrid1D`](@ref).
"""
function SectorWilson1D(H, NX::Int, O::AbstractMatrix, bandindices::AbstractArray; kwargs...)
    isempty(bandindices) && throw(ArgumentError("select the occupied bandindices explicitly"))
    NX > 1 || throw(ArgumentError("a Wilson grid needs at least two points"))
    _, statesgrid = statesgrid1D(H, NX, bandindices; kwargs...)
    sectorwilson(statesgrid, O)
end

"""
    statesgrid1Dslice(H, NX::Int, bandindices, dim::Int, kfixed::Real; kwargs...)

Eigenvectors along reduced momentum direction `dim` of a 2D Hamiltonian,
with the other coordinate fixed at `kfixed`. For `dim=1`, sample `[x, kfixed]`;
for `dim=2`, sample `[kfixed, x]`. Delegates to [`statesgrid1D`](@ref), including
Hamiltonian validation and execution keywords, and returns `(kgrid, states)`.
"""
function statesgrid1Dslice(H, NX::Int, bandindices::AbstractArray, dim::Int, kfixed::Real; kwargs...)
    dim in (1, 2) || throw(ArgumentError("dim must be 1 or 2"))
    Hslice = k -> H(dim == 1 ? [k[1], kfixed] : [kfixed, k[1]])
    statesgrid1D(Hslice, NX, bandindices; kwargs...)
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

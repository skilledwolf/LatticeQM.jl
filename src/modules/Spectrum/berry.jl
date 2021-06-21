import LinearAlgebra: det, norm, I, svd
import SharedArrays: SharedArray

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

    for i=1:n-1, j=1:n-1
        S00 = statesgrid[i,  j, :, :]
        S10 = statesgrid[i+1,j, :, :]
        S01 = statesgrid[i, j+1, :, :]
        S11 = statesgrid[i+1, j+1, :, :]
        F[i,j] = plaquettephase(S00, S10, S01, S11)
    end

    F
end

"""
    statesgrid(H, NX::Int, NY::Int=0, bandindices::AbstractArray=[])

Evaluates the eigenvectors on a discretized grid (2D Hamiltonian only!) and stores the result (preserving the grid information).
This method is useful when plaquette phases need to be calculated.
"""
function statesgrid(H, NX::Int, NY::Int=0, bandindices::AbstractArray=[])
    # Prepare indices and sizes
    NY = (NY<1) ? NX : NY
    indices = collect(Iterators.product(1:NX, 1:NY))

    function wavefunctionsf(k)
        wavefunctions(H(k))
    end
    
    M1 = size(wavefunctionsf(zeros(2)), 2) # dimension of Hilbert space
    bandindices = (bandindices==[]) ? collect(1:M1) : bandindices
    M2 = size(bandindices,1) # number of occupied bands

    # Prepare k-grid
    kgrid = [[x;y] for x=range(0; stop=1, length=NX), y=range(0; stop=1, length=NY)]
    midkgrid = [[1/(2*NX)+x;1/(2*NY)+y] for x=range(0; stop=1-1.0/NX, length=NX-1), y=range(0; stop=1-1.0/NY, length=NY-1)]

    # Compute eigenspectrum on the k-grid
    statesgrid0 = convert(SharedArray, zeros(ComplexF64, NX, NY, M1, M2))
    @sync @distributed for (i_,j_) in indices
        statesgrid0[i_,j_, :, :] = wavefunctionsf(kgrid[i_,j_].+1.34e-8)[:,bandindices]
    end

    kgrid, midkgrid, statesgrid0
end

"""
    berry(H, NX::Int, NY::Int=0, bandindices::AbstractArray=[])

Convenience method for `berry(statesgrid)`.
"""
function berry(H, NX::Int, NY::Int=0, bandindices::AbstractArray=[1])
    # wavefunctions(k::Vector) -> Matrix{Complex}
    kgrid, midkgrid, statesgrid0 = statesgrid(H, NX, NY, bandindices)
    midkgrid, berry(statesgrid0)
end

import ..Structure: rotation2D

"""
    berryalongpath(H, kpoints)

Calculate the abelian Berry curvature for each band along a path of discrete k points.
It builds little plaquettes along the path between the kpoints[:,i] and kpoints[:,i+1].
"""
function berryalongpath(H, kpoints)

    function wavefunctionsf(k)
        wavefunctions(H(k))
    end

    N = size(kpoints,2) # number of k points
    M = size(wavefunctionsf(zero(first(eachcol(kpoints)))), 2) # number of bands

    # add some "dummy" points to build the grid
    k0 = kpoints[:,1] - (kpoints[:,2]-kpoints[:,1])
    kNp1 = kpoints[:,N] + (kpoints[:,N]-kpoints[:,N-1])
    kpoints = hcat(k0,kpoints[:,:],kNp1)

    kpoints .+= 0.10*norm(kpoints[:,2]-kpoints[:,1])*rand(2)

    # prepare the array for the results
    berryc = convert(SharedArray, zeros(M, N))

    # @sync @distributed for j_=2:N+1
    pmap(2:N+1) do j_

        k0 = kpoints[:,j_]

        δkR = (kpoints[:,j_+1]-k0)/2
        δkL = (kpoints[:,j_-1]-k0)/2
        δkU = (rotation2D(π/2)*δkR + rotation2D(-π/2)*δkL)/2
        δkD = (rotation2D(-π/2)*δkR + rotation2D(π/2)*δkL)/2

        plaquette = k0 .+  hcat(δkR, δkU, δkD, δkL)

        States = [wavefunctionsf(k) for k in eachcol(plaquette)]

        for i_=1:M
            berryc[i_,j_-1] = plaquettephase(
                (U[:,i_] for U in States)...
            )
        end

        nothing
    end

    berryc
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


function statesgrid1D(H, NX::Int, bandindices::AbstractArray=[])
    function wavefunctionsf(k)
        wavefunctions(H(k))
    end

    M1 = size(wavefunctionsf(zeros(1)), 2) # dimension of hilbert space
    bandindices = (bandindices==[]) ? collect(1:M1) : bandindices
    M2 = size(bandindices,1) # number of occupied bands

    # Prepare k-grid
    kgrid = LinRange(0,1,NX)

    # Compute eigenspectrum on the k-grid
    statesgrid0 = convert(SharedArray, zeros(ComplexF64, NX, M1, M2))
    @sync @distributed for i_ in 1:NX
        statesgrid0[i_, :, :] = wavefunctionsf([kgrid[i_]])[:,bandindices]
    end

    kgrid, statesgrid0
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

    en, U = spectrum(F)
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
    cherns = [sum(Spectrum.berry(grid[:,:,:,i:i])[2]) for i=1:size(grid,4)]

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

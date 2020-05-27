function L(State1::T, State2::T) where {T<:AbstractArray{<:Complex,N}} where N
    res = det(State1' * State2)
    return res/abs(res)
end

@legacyalias plaquettephase PlaquettePhase
function plaquettephase(S00::T, S10::T, S01::T, S11::T) where {T<:AbstractArray{<:Complex, N}} where N
    real(1.0/(2π*1.0im) * log( L(S00, S10) * L(S10, S11) * L(S01, S11)^(-1) * L(S00,S01)^(-1) ))
end

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

@legacyalias berry BerryF
function berry(wavefunctions::Function, NX::Int, NY::Int=0, bandindex=[1])
    # wavefunctions(k::Vector) -> Matrix{Complex}

    if NY < 1
        NY = NX
    end

    M1 = size(wavefunctions(zeros(2)), 2)
    M2 = size(bandindex,1)

    kgrid = [[x;y] for x=range(0; stop=1, length=NX), y=range(0; stop=1, length=NY)]

    statesgrid = convert(SharedArray, zeros(ComplexF64, NX, NY, M1, M2))

    indices = collect(Iterators.product(1:NX, 1:NY))

    @sync @distributed for (i_,j_) in indices
        statesgrid[i_,j_, :, :] = wavefunctions(kgrid[i_,j_])[:,bandindex]
    end

    berry(statesgrid)
end

@legacyalias berryalongpath BerryCurvature
function berryalongpath(wavefunctions::Function, kpoints::AbstractMatrix{<:Float64})
"""
    Calculate the abelian Berry Curvature for each band along a path of discrete k points
"""

    N = size(kpoints,2) # number of k points
    M = size(wavefunctions(zero(first(eachcol(kpoints)))), 2) # number of bands

    # add some "dummy" points to build the grid
    k0 = kpoints[:,1] - (kpoints[:,2]-kpoints[:,1])
    kNp1 = kpoints[:,N] + (kpoints[:,N]-kpoints[:,N-1])
    kpoints = hcat(k0,kpoints,kNp1)

    kpoints .+= 0.10*norm(kpoints[:,2]-kpoints[:,1])*rand(2)

    # prepare the array for the results
    berryc = convert(SharedArray, zeros(M, N))

    @sync @distributed for j_=2:N+1

        k0 = kpoints[:,j_]

        δkR = (kpoints[:,j_+1]-k0)/2
        δkL = (kpoints[:,j_-1]-k0)/2
        δkU = (rotation2D(π/2)*δkR + rotation2D(-π/2)*δkL)/2
        δkD = (rotation2D(-π/2)*δkR + rotation2D(π/2)*δkL)/2

        plaquette = k0 .+  hcat(δkR, δkU, δkD, δkL)

        States = [wavefunctions(k) for k in eachcol(plaquette)]

        for i_=1:M
            berryc[i_,j_-1] = plaquettephase(
                [U[:,i_] for U in States]...
            )
        end

    end

    berryc
end

@legacyalias getberry! get_berry!
function getberry!(bands::BandData, h, ks)
    obs = Array(berryalongpath(wavefunctions(h), ks.points))
    obs = reshape(obs, (size(obs)...,1))

    if bands.obs == nothing
        bands.obs = obs
    else
        bands.obs = cat(bands.obs, obs, dims=(3))
    end

    nothing
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

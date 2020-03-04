using SparseArrays

function graphene(lat::Lattice; mode=:nospin, format=:auto, cellrange=2, kwargs...)

    t(args...) = t_graphene(args...; kwargs...)

    hops = gethops(lat, t; cellrange=cellrange, format=format, vectorized=true)

    if mode==:spinhalf
        hops = addspin(hops, mode)
    end

    hops
end


valleyoperator(args...; kwargs...) = valleyoperator!(Hops(), args...; kwargs...)
function valleyoperator!(hops, lat::Lattice; spinhalf=false, zmode=:anti, kwargs...)
    @assert latticedim(lat) == 2 && countorbitals(lat) > 1
    addhaldane!(hops, lat, √3/9; spinhalf=spinhalf, mode=:anti, zmode=zmode, kwargs...)
end

@legacyalias addsublatticeimbalance! add_sublatticeimbalance
function addsublatticeimbalance!(hops, lat::Lattice, Δ::AbstractFloat; kwargs...)

    # Only go through the trouble of constructing this matrix for finite Δ
    if abs(Δ) ≈ 0
        return nothing
    end

    μ = Δ .* (extrapositions(lat, "sublattice") .- 0.5)
    addchemicalpotential!(hops, lat, vec(μ))

    nothing
end

@legacyalias addhaldane add_haldane
function addhaldane!(hops, lat::Lattice, t2; ϕ=π/2, spinhalf=false, mode=:none, zmode=:none)
    """
    This method is a somewhat inefficient way to compute the haldane hopping matrix.
    The only upside to it is that it uses methods that I already implemented and
    that it is fairly general.
    """

    t2 = hcat(t2) # turn scalars into 1x1 matrix

    if norm(t2) ≈ 0
        return nothing
    end

    if spinhalf
        t2 = addspin(t2, :spinhalf)
    end

    d = size(t2,1)

    cross2D(vec1, vec2) = vec1[1] * vec2[2] - vec1[2] * vec2[1] # needed later on in this scope

    # NN  = find_neighbors(lat, 1.0)
    NNN = Structure.getneighbors(lat, √3)

    N = countorbitals(lat)
    R = positions(lat) # positions of atoms within unit cell
    sublattice = extrapositions(lat, "sublattice") .- 0.5
    if hasdimension(lat, "z")
        zpositions = extrapositions(lat, "z")
    end
    A = getA(lat)

    neighbors = [[i;j] for i=-1:1 for j=-1:1]
    δAs = [A * v for v in neighbors]

    # hops = Dict{Vector{Int},SparseMatrixCSC{ComplexF64}}()
    for δR = keys(NNN)
        if !haskey(hops, δR)
            hops[δR] = spzeros(ComplexF64, N*d, N*d)
        end
    end

    for (δR, NNN_pairs) = NNN
        for (i,j) in NNN_pairs

            Ri = R[:,i] .+ (A * δR)
            Rj = R[:,j]

            if mode==:sublatticeA
                s = (sublattice[i]>0) ? 1 : 0 #sign(sublattice[i])
            elseif mode==:sublatticeB
                s = (sublattice[i]<0) ? 1 : 0
            elseif mode==:anti
                s = sign(sublattice[i])
            else
                s = 1
            end

            if hasdimension(lat, "z")
                if zmode==:positive
                    s *= (zpositions[i]>0) ? 1 : 0 #(zpositions[i] > 0) ? 1 : 0
                elseif zmode==:negative
                    s *= (zpositions[i]<0) ? 1 : 0
                elseif zmode==:anti
                    s *= sign(zpositions[i]+1e-4)
                end
            end

            for δAi=δAs, δri=eachcol(R)

                R0 = δAi .+ δri

                if 0.9 < norm(Ri-R0) < 1.1 && 0.9 < norm(Rj-R0) < 1.1 # we found the common
                    hops[δR][1+d*(i-1):d+d*(i-1), 1+d*(j-1):d+d*(j-1)] += s * t2 * I * 1.0im * sign( cross2D(R0.-Rj, Ri.-R0) ) #* exp(1.0im * ϕ * sign( cross2D(R0.-Rj, Ri.-R0) ) )
                    break
                end
            end
        end
    end

    hops
end

@legacyalias gethaldane get_haldane_hops
function gethaldane(args...; kwargs...)
    newhops = Hops()
    addhaldane!(newhops, args...; kwargs...)
end

@legacyalias addspinorbit! add_spinorbit!
function addspinorbit!(hops, lat, t2, args...)

    if t2≈0.0
        return nothing
    end

    newhops = Hops()
    addhaldane!(newhops, lat, t2, args...)
    addhops!(hops, kron(newhops, σZ))
end

using ..TightBinding: MAX_DENSE, MAX_DIAGS

@legacyalias getrashba rashba_hops
function getrashba(lat::Lattice, λ::Function; cellrange=2, format=:auto, tmin=1e-7)

    function hop(R1, R2)
        # this is a "vectorized" hopping function, i.e., it expects
        # coordinate matrices R1 and R2 as inputs
        # (had to be implemented this way for performance)

        N = size(R1,2)
        d = 2

        # Preallocate memory: important for huge sparse matrices
        maxind = (N>MAX_DENSE) ? round(Int, MAX_DIAGS * N) : N^2 # MIN_SPARSITY * N^2 # semi-arbitrary limit for dense allocation
        IS = Vector{Int}(undef, maxind*d^2)
        JS = similar(IS)
        VS = similar(IS, ComplexF64)
        δR = similar(R1)
        R0 = similar(R1)

        count = 0
        # do not tamper with this loop.
        # this way is orders of magnitude faster than any other implementation
        @fastmath @inbounds for j=1:N
            @views δR .= R1 .- R2[:,j]
            @views R0 .= (R1 .+ R2[:,j])./2

            for i=1:N
                @views Δ = sqrt(sum(abs2,δR[1:3,i]))

                if 1.1 < Δ || Δ < 0.9 || abs(δR[3,i]) > 0.2
                    continue
                end

                V = 1.0im .* λ(R0[:,i]) .* (σX.*δR[2,i] .- σY.*δR[1,i])

                for i0=1:d, j0=1:d
                    if abs(V[i0,j0]) < tmin
                        continue
                    end
                    count = count+1
                    IS[count], JS[count], VS[count] = (i-1)*d + i0,  (j-1)*d + j0, V[i0,j0]
                end
            end
        end

        @views sparse(IS[1:count],JS[1:count],complex(VS[1:count]), d*N, d*N)
    end

    gethops(lat, hop; vectorized=true, cellrange=cellrange, format=format)
end
getrashba(lat::Lattice, λ::AbstractFloat; kwargs...) = getrashba(lat, x->λ; kwargs...)

@legacyalias addrashba! add_rashba!
addrashba!(hops, lat, rashba::Function; kwargs...) = addhops!(hops, getrashba(lat, rashba; kwargs...))
function addrashba!(hops, lat, rashba::Number)
    if !(rashba≈0.0)
        addrashba!(hops, lat, x->rasbha)
    end
end

################################################################################
################################################################################
################################################################################

"""
    This is the hopping amplitude between to twisted layers as given by overlap
    intergrals between p_z orbitals.
"""

norm2(r) = dot(r,r)

@legacyalias t_graphene graphene_hops

using ..Utils: heaviside

graphene_interlayer(δz, Δ; z, ℓinter) =  δz^2 /Δ^2 * exp(-(Δ-z)/ℓinter)
graphene_intralayer(δz, Δ; a, ℓintra, ℓz) = (1 - δz^2 /(Δ^2)) * exp(-(Δ-a)/ℓintra -δz^2 /ℓz^2)

@inline function t_graphene(r::AbstractVector{Float64}; tz::Float64=0.46, t0::Float64=1.0, ℓinter::Float64=0.125, ℓintra::Float64=0.08, ℓz::Float64=0.01,z::Float64=3.0, a::Float64=1.0,
           Δmin::Float64=0.1, Δmax::Float64=5.0)
    @views Δ = sqrt.(sum(abs2,r[1:3]))
    result = 0.0
    if Δmax > Δ > Δmin
        @views δz = r[3]
        χ = δz^2 /(Δ^2)
        result +=  (-t0) * (1-χ) * exp(-(Δ-a)/ℓintra) * exp(-δz^2 /ℓz^2) - tz * χ * exp(-(Δ-z)/ℓinter)
    end
    result
end

using ..TightBinding: MAX_DENSE, MAX_DIAGS

function t_graphene(R1::Matrix{Float64}, R2::Matrix{Float64}; tmin=1e-7, tz::Float64=0.46, t0::Float64=1.0,
    ℓinter::Float64=0.125, ℓintra::Float64=0.08, ℓz::Float64=0.001,z::Float64=3.0, a::Float64=1.0,
    Δmin::Float64=0.1, Δmax::Float64=5.0,
    kwargs...)

    N = size(R1,2)

    # Preallocate memory: important for huge sparse matrices
    maxind = (N>MAX_DENSE) ? round(Int, MAX_DIAGS * N) : N^2 # MIN_SPARSITY * N^2 # semi-arbitrary limit for dense allocation
    IS = Vector{Int}(undef, maxind)
    JS = similar(IS)
    VS = similar(IS, Float64)
    δR = similar(R1)

    count = 0
    @fastmath @inbounds for j=1:N
        @views δR .= R1 .- R2[:,j]

        for i=1:N
            @views Δ = sqrt(sum(abs2,δR[1:3,i])) # absolute value of distance vector

            if Δmax < Δ || Δ < Δmin
                continue
            end

            δz = δR[3,i]
            χ = δz^2 /(Δ^2)

            if δz < a
                v = -t0 * (1-χ) * exp(-(Δ-a)/ℓintra) #* exp(-δz^2 /ℓz^2)
            else
                v = -tz * χ * exp(-(Δ-z)/ℓinter)
            end

#             v =  -t0 * (1-χ) * exp(-(Δ-a)/ℓintra) * exp(-δz^2 /ℓz^2) - tz * χ * exp(-(Δ-z)/ℓinter)

            if abs(v) < tmin
                continue
            end

            count = count+1
            IS[count], JS[count], VS[count] = i, j, v
        end
    end

    @views sparse(IS[1:count],JS[1:count],complex(VS[1:count]), N, N)
end

# ####
# # This version is more natural but 4x slower than the fastest implementation
# ####
# function t_graphene(R1::Matrix{Float64}, R2::Matrix{Float64}; tmin=1e-7,
#            kwargs...) #where {T<:AbstractMatrix}
#
#     N = size(R1,2)
#
#     # Preallocate memory: important for huge sparse matrices
#     maxind = (N>MAX_DENSE) ? round(Int, MAX_DIAGS * N) : N^2 # semi-arbitrary limit for dense allocation
#     IS = Vector{Int}(undef, maxind)
#     JS = similar(IS)
#     VS = similar(IS, Float64)
#
#     δR = similar(R1)
#
#     t0(args...) = t_graphene(args...; kwargs...)
#
#     count = 1
#     @fastmath @inbounds for j=1:N
#         @views δR .= R1 .- R2[:,j]
#         for i=1:N
#             @views v = t0(δR[:,i]) #1.4e-7 * rand() # t(δr)
#             if abs(v) > tmin
#                 IS[count], JS[count], VS[count] = i, j, v
#                 count = count+1
#             end
#         end
#     end
#     count = count - 1
#
#
#     @views sparse(IS[1:count],JS[1:count],complex(VS[1:count]), N, N)
# end
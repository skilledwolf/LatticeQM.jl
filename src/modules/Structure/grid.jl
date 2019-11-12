
using RecursiveArrayTools

function regulargrid(;nk::Int=100, dim::Int=2)
    nk = floor(Int, sqrt(nk))

    it1d = range(0.0, 1.0; length=nk+1)[1:end-1]
    itNd = Iterators.product(Iterators.repeated(it1d, dim)...)

    out = convert(Array, VectorOfArray([[y...] for y=[itNd...]]))

    out
end

rot(θ::Float64) = [cos(θ) -sin(θ); sin(θ) cos(θ)]

function randomgrid(;nk::Int=100, dim::Int=2, rot_symmetry::Int=1, B=:id)

    if B==:id
        B = Matrix{Float64}(I, dim, dim)
    end

    @assert rot_symmetry>0

    N = div(nk, rot_symmetry) # integer division
    @assert N>0

    @info("Random k grid with (symmetrized) points", N*rot_symmetry)
    ks = rand(Float64, (dim,N))

    # symmetrized sampling
    for i=1:(rot_symmetry-1)
        ks = Matrix([ks inv(B)*rot(2π/rot_symmetry*i)*B*ks])
    end

    ks
end

function foldBZ!(M, kpoints::AbstractMatrix)
    """
    M = B^T B, where B=[G1 G2] contains the reciprocal lattice vectors as columns
    kpoints has the lattice points as columns in units of lattice vectors
    """
    G0sq = sum(M)
    b = M[1,2]/G0sq

    for j_ = 1:size(kpoints,2)
        kpoints[:,j_] .= mod.(kpoints[:,j_], 1.0)
        α = M * kpoints[:,j_]

        α1 = α[1]/M[1,1]
        α2 = α[2]/M[2,2]
        α3 = (α[1]+α[2])/G0sq

        if α1 > 1/2 && α2 < 1/2+b
            δk = [1.0, 0.0]
        elseif α1 < 1/2+b && α2 > 1/2
            δk = [0.0, 1.0]
        elseif α3 > 1/2
            δk = [1.0, 1.0]
        else
            δk = [0.0,0.0]
        end
        kpoints[:,j_] -= δk
    end

    nothing
end
foldBZ!(lat::Lattice, kpoints::AbstractMatrix) = foldBZ!(transpose(get_B(lat))*get_B(lat),kpoints)
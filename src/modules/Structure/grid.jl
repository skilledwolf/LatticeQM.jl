
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

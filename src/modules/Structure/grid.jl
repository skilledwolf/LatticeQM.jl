
function regulargrid(;nk::Int=100, dim::Int=2)
    nk = floor(Int, nk^(1/dim))

    it1d = range(0.0, 1.0; length=nk+1)[1:end-1]
    itNd = Iterators.product(Iterators.repeated(it1d, dim)...)

    out = hcat(([v...] for v=itNd)...)

    out
end

function randomgrid(;nk::Int=100, dim::Int=2, rot_symmetry::Int=1, B=:id)
    @assert rot_symmetry>0

    if B==:id
        B = Matrix(1.0*I, dim, dim)
    end

    N = div(nk, rot_symmetry) # integer division
    @assert N>0

    println("Random (symmetrized) k grid. # points: ", N*rot_symmetry)
    ks = rand(Float64, (dim,N))

    # symmetrized sampling
    for i=1:(rot_symmetry-1)
        ks = Matrix([ks inv(B[1:dim,1:dim])*rotation2D(2π/rot_symmetry*i)*B[1:dim,1:dim]*ks])
    end

    ks
end
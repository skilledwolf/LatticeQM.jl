

function regulargrid(N::Int, dim::Int=2)
    it1d = range(0.0, 1.0; length=N)
    itNd = Iterators.product(Iterators.repeated(it1d, dim)...)

    out = convert(Array, VectorOfArray([[y...] for y=[itNd...]]))

    out
end

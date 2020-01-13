using ..TightBinding

function initialguess(v::AnyHops, mode=:random; lat=:nothing)
    N = size(first(values(v)), 1)

    ρs = Hops()
    for δL=keys(v)
        ρs[δL] = zeros(ComplexF64, size(v[δL]))
    end

    if mode==:randombig
        mat = rand(ComplexF64, N, N)
        ρs[zero(first(keys(ρs)))] = (mat + mat') ./ 2

    elseif mode==:random
        @assert mod(N,2)==0
        n = div(N,2)
        mat = rand(ComplexF64, n, n)

        # Generate a random spin orientation at a lattice site
        function randmat()
            # d = rand(Float64, 3)
            d = -1.0 .+ 2 .* rand(Float64, 3)
            p = 0.5 .* (σ0 .+ sum(d[i_]/norm(d) .* σs[i_] for i_=1:3))
        end

        ρs[zero(first(keys(ρs)))] = Matrix(sum(kron(sparse([i_],[i_], [1.0+0.0im], n,n), randmat()) for i_=1:n))

    elseif mode==:antiferro || mode==:antiferroZ
        sublA, sublB = getoperator(lat, ["sublatticeA", "sublatticeB"])
        mat = sublA .- sublB

        σUP = 0.5 .* (σ0 .+ σZ)

        ρs[zero(first(keys(ρs)))] = kron(mat,σUP)

    elseif mode==:antiferroX
        sublA, sublB = getoperator(lat, ["sublatticeA", "sublatticeB"])
        mat = sublA .- sublB

        σUP = 0.5 .* (σ0 .+ σX)

        ρs[zero(first(keys(ρs)))] = kron(mat,σUP)

    elseif mode==:ferro || mode==:ferroZ #|| mode==:ferroz
        @assert mod(N,2)==0
        n = div(N,2)
        σUP = 0.5 .* (σ0 .+ σZ)
        ρs[zero(first(keys(ρs)))] =  2. * kron(Diagonal(ones(n)),σUP)

    elseif mode==:ferroX #|| mode==:ferroz
        @assert mod(N,2)==0
        n = div(N,2)
        σLEFT = 0.5 .* (σ0 .+ σX)
        ρs[zero(first(keys(ρs)))] =  2. * kron(Diagonal(ones(n)), σLEFT)

    else
        error("Unrecognized mode '$mode' in initialguess(...).")

    end

    ρs
end


###################################################################################################
# Backwards compatibility
###################################################################################################
export initial_guess
@legacyalias initialguess initial_guess

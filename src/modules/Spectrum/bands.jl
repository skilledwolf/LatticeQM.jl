using Distributed
using ProgressMeter

function bandmatrix(h::Function, ks::AbstractMatrix; num_bands=nothing, kwargs...)
    L = size(ks)[2]         # no. of k points

    if num_bands != nothing
        M = num_bands
    else
        M = size(h(ks[:,1]))[1] # no. of bands
    end


    bands = convert(SharedArray, zeros(Float64, M, L))

    energies = ϵs(h; num_bands=num_bands, kwargs...)
    @sync @showprogress 1 "Computing bands..."  @distributed for j_=1:L

        bands[:,j_] .= real.(energies(ks[:,j_]))
    end

    Matrix(bands)
end

coloredbandmatrix(h::Function, ks::AbstractMatrix, projector::Function; kwargs...) = coloredbandmatrix(h,ks,[projector];kwargs...)
function coloredbandmatrix(h::Function, ks0::AbstractMatrix, projector::AbstractVector; num_bands=nothing, kwargs...)
    """
        h(k): returns hermitian Hamiltonian at k-point
        ks: collection of k points (see kIterable)
        projector: function that returns a real value for a given (k,ψ,E_k).
            can also be as Vector of such functions.
    """
    N = size(ks0, 2)        # no. of k points
    M = (num_bands!=nothing) ? num_bands : size(h(ks0[:,1]), 1)  # no. of bands

    bands = zeros(Float64, M, N)
    obs = zeros(Float64, M, N, length(projector))

    Σ = spectrum(h; num_bands=num_bands, kwargs...)
    projector = projector

    # Parallized loop
    bands = convert(SharedArray, bands)
    obs = convert(SharedArray, obs)

    p = Progress(N, "Computing bands... ")
    channel = RemoteChannel(()->Channel{Bool}(N), 1)

    @sync begin
        # this task prints the progress bar
        @async begin
            done = 0
            while done < N
                take!(channel)
                next!(p)
                done = done + 1
            end
        end

        # this task does the computation
        @async begin
            @distributed for j_=1:N   #@sync @distributed for j_=1:N
                k = ks0[:,j_]
                ϵs, U = Σ(k)

                bands[:,j_] .= real.(ϵs)

                for (i_,ψ)=enumerate(eachcol(U))
                    for (n_, proj)=enumerate(projector)
                        obs[i_,j_,n_] = proj(k,ψ,ϵs[i_])
                    end
                end

                put!(channel, true)
            end
#             put!(channel, false) # this tells the printing task to finish
        end
    end

    bands = convert(Array, bands)
    obs = convert(Array, obs)

    bands, obs
end

function get_bands(h,ks; projector=nothing, kwargs...)
    if projector != nothing
        return coloredbandmatrix(h,ks,projector; kwargs...)
    else
        bands = bandmatrix(h,ks; kwargs...)
        return bands, nothing
    end
end

function get_bands(h::Function, ks::DiscretePath; kwargs...)

    bands, obs = get_bands(h, points(ks); kwargs...)

    BandData(bands, obs, ks)
end

"""
    bandmatrix(H, ks, drive, M::Integer; num_bands::Int=0, kwargs...)
Returns all(!) Floquet bands of a system with time-averaged Hamiltonian H, driven by a periodic potential given in drive.
The path is chosen along ks and M denotes the cutoff frequency for not regarding Fourier modes anymore.
A d-dimensional system with a cutoff M will lead to d*(2*M+1) bands.
"""
function Spectrum.bandmatrix(H, ks, drive, M::Integer; kwargs...)

    HFloq = getFloquetHamiltonian(M, Spectrum.ensureH(H), drive) #new
    return Spectrum.bandmatrix(HFloq, ks; kwargs...)
end

# """Same as bandmatrix but allows for multithreading."""
# function Spectrum.bandmatrix_multithread(H, ks, drive::periodicDrive, M::Integer; num_bands::Int=0, kwargs...)
    
#     HFloq = getFloquetHamiltonian(M, Spectrum.ensureH(H), drive)
#     return Spectrum.bandmatrix_multithread(HFloq, ks; kwargs...)
# end

"""Version of bandmatrix with observables. Returns stroboscopic expectation values! (not the same as time averaged expectation values!)"""
function Spectrum.bandmatrix(H, ks, projector, drive::periodicDrive, M::Integer; num_bands::Int=0, kwargs...)
    projector = Spectrum.handleprojector(projector)
    ks = points(ks)
    dimH = dim(H, ks) #dimension of the Bloch Hamiltonian
    HFloq = getFloquetHamiltonian(M, Spectrum.ensureH(H), drive)
    if !(num_bands>0)
        num_bands = dim(HFloq, ks) #number of total bands, in the end we reduce to the frist Floquet BZ
    end

    N = size(ks, 2) # number of k points
    L = length(projector)
    #bands = convert(SharedArray, zeros(Float64, num_bands, N))
    #obs   = convert(SharedArray, zeros(Float64, num_bands, N, L))
    bands = zeros(Float64, num_bands, N)
    obs = zeros(Float64, dimH, N, L)

    spectrumf = spectrum(HFloq; num_bands=num_bands, format=:sparse, kwargs...)

    @sync @showprogress 1 "Computing bands... " @distributed for j_=1:N
#     @showprogress 1 "Computing bands..." for j_=1:N
        ϵs, U = spectrumf(ks[:,j_])
        bands[:,j_] .= real.(ϵs)
        U = transformeigvecs(U, M, 0) #evaluate at t=0

        for i_=1:dimH, n_=1:L
            obs[i_,j_,n_] = projector[n_](ks[:,j_],U[:,i_],ϵs[i_])
        end
    end

    Array(bands), Array(obs)
end

#### Floquet specific utilities ###############################################################################
"""function reducetoFBZ(bands::AbstractArray, M::Integer)
    Takes a two dimensional array denoting a list of bands which came from diagonalizing a Floquet Hamiltonian
    and returns the bands lying in the first Floquet Brioullin zone."""
function reducetoFBZ(bands::AbstractArray, M::Integer)
    dim = Int(size(bands,1)/(2*M+1)) #dimension of the Bloch Hamiltonian
    Fbands = zeros(dim,size(bands,2))

    for i=1:dim
        Fbands[i,:] = bands[M*dim+i,:]
    end

    return Fbands
end

"""transformeigvecs(U::AbstractMatrix, M::Integer, t::Number)
    Converts a matrix of eigenvectors U obtained from solving the infinite (truncated at M)
    Floquet problem to eigenvectors in time space (at time t)."""
function transformeigvecs(U::AbstractMatrix, M::Integer, t::Number)
    dim = Int(size(U,1)/(2*M+1)) #dimension of the Bloch Hamiltonian
    Unew = zeros(ComplexF64,dim, size(U,2))

    for i=1:size(U,2)
        Unew[:,i] = transform(U[:,i], dim, M, t)
    end

    return Unew
end

"""Converts an eigenvector vec obtained from solving the infinite (truncated at M)
    Floquet problem to an eigenvector in time space (at time t)."""
function transform(vec::Array{ComplexF64,1}, dim::Integer, M::Integer, t::Number)
    vecnew = zeros(ComplexF64,dim)

    for i=0:(2*M)
        vecnew += exp(-im*(i-M)*t).*vec[i*dim+1:(i+1)*dim]
    end

    return vecnew
end


"""From bands.jl in LatticeQM: I added additional functions taking a periodicDrive as argument.
Maybe it would be wise to put it into that file in the end."""
function Spectrum.getbands_parallel(H, ks::DiscretePath, drive, M::Integer;kwargs...)
    bands = Spectrum.bandmatrix(H, points(ks), drive, M; kwargs...)
    bands = reducetoFBZ(bands, M) # reduce the bands to the first Floquet-Brioullin-zone. We choose the result from the middle of the spectrum because they are the most accurate.
    obs = nothing
    Spectrum.BandData(bands, obs, ks)
end

function Spectrum.getbands_parallel(H, ks::DiscretePath, projector, drive::periodicDrive, M::Integer; kwargs...)
    bands, obs = Spectrum.bandmatrix(H, points(ks), projector, drive, M; kwargs...)
    bands = reducetoFBZ(bands, M) # reduce the bands to the first Floquet-Brioullin-zone.
    Spectrum.BandData(bands, obs, ks)
end


# function Spectrum.getbands_multithread(H, ks::DiscretePath, drive, M::Integer; kwargs...)
#     bands = Spectrum.bandmatrix_multithread(H, points(ks), drive, M; kwargs...)
#     bands = reducetoFBZ(bands, M) # reduce the bands to the first Floquet-Brioullin-zone. We choose the result from the middle of the spectrum because they are the most accurate.
#     obs = nothing
#     Spectrum.BandData(bands, obs, ks)
# end

# function Spectrum.getbands_multithread(H, ks::DiscretePath, projector, drive::periodicDrive, M::Integer; kwargs...)
#     bands, obs = Spectrum.bandmatrix_multithread(H, points(ks), projector, drive, M; kwargs...)
#     bands = reducetoFBZ(bands, M) # reduce the bands to the first Floquet-Brioullin-zone.
#     Spectrum.BandData(bands, obs, ks)
# end


"""
    wavefunctions(H, drive::periodicDrive, M::Integer; timespace=false)

Returns a matrix of Floquet eigenstates of a system with time-averaged Hamiltonian H,
driven by the periodic potential contained in drive.
M denotes the frequency cutoff and therefore also the length of the eigenvectors in frequency space.
For a d-dimensional Hamiltonian, the Fourier-space eigenvectors will have length d*(2*M+1).
If timespace=true, these eigenvectors are converted back to time space and will have length dim.
"""
function Spectrum.wavefunctions(H::AbstractMatrix, drive::periodicDrive, M::Integer; timespace=false)
    HFloq = getFloquetHamiltonian(M, Spectrum.ensureH(H), drive)
    U = Spectrum.wavefunctions(HFloq)

    if timespace
        U = transformeigvecs(U, M, 0)
    end

    return U
end

Spectrum.wavefunctions(H::Function, drive::periodicDrive, M::Integer; kwargs...) = k -> Spectrum.wavefunctions(H(k), drive::periodicDrive, M::Integer; kwargs...)
Spectrum.wavefunctions(H::Function, drive::Function, M::Integer; kwargs...) = k -> Spectrum.wavefunctions(H(k), drive(k), M::Integer; kwargs...)

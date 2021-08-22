

#### Floquet specific utilities ###############################################################################
"""function reducetoFBZ(bands::AbstractArray, M::Integer)
    Takes a two dimensional array denoting a list of bands which came from diagonalizing a Floquet Hamiltonian
    and returns the bands lying in the first Floquet Brioullin zone."""
function getfirstFBZ(bands::AbstractArray, M::Integer)
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


import ..Structure.Paths

"""
    keepfirstFBZ!(data::Spectrum.BandData, H::FloquetOperator)

Throw away all band data except for the first Floquet Brillouin zone.
"""
function keepfirstFBZ!(data::Spectrum.BandData, H::FloquetOperator)
    bands = getfirstFBZ(data.bands, H.M)
    data.bands = bands
    data
end


# """
#     wavefunctions(H, drive::periodicDrive, M::Integer; timespace=false)

# Returns a matrix of Floquet eigenstates of a system with time-averaged Hamiltonian H,
# driven by the periodic potential contained in drive.
# M denotes the frequency cutoff and therefore also the length of the eigenvectors in frequency space.
# For a d-dimensional Hamiltonian, the Fourier-space eigenvectors will have length d*(2*M+1).
# If timespace=true, these eigenvectors are converted back to time space and will have length dim.
# """
# function Spectrum.wavefunctions(HF::FloquetOperator; timespace=false)

#     HFloq = getFloquetHamiltonian(HF)
#     U = Spectrum.wavefunctions(HFloq)

#     if timespace
#         U = transformeigvecs(U, HF.M, 0)
#     end

#     return U
# end



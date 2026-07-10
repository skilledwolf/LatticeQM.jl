

#### Floquet specific utilities ###############################################################################
"""getfirstFBZ(bands::AbstractArray, M::Integer, omega::Real)
    Takes a two dimensional array denoting a list of bands which came from diagonalizing a Floquet Hamiltonian
    and returns the bands lying in the first Floquet Brillouin zone, i.e., for each k point the dim
    quasienergies inside the window (-omega/2, omega/2], where dim is the dimension of the Bloch Hamiltonian.
    If the number of eigenvalues inside the window differs from dim (degeneracies at the zone edge or
    truncation artifacts), the dim states closest to the window are selected instead and folded back
    into (-omega/2, omega/2]."""
function getfirstFBZ(bands::AbstractArray, M::Integer, omega::Real)
    dim = Int(size(bands,1)/(2*M+1)) #dimension of the Bloch Hamiltonian
    Fbands = zeros(dim,size(bands,2))

    halfomega = omega/2
    inwindow(e) = -halfomega < e <= halfomega
    fold(e) = inwindow(e) ? e : halfomega - mod(halfomega - e, omega) #fold into (-omega/2, omega/2]
    windowdist(e) = e > halfomega ? e - halfomega : (e <= -halfomega ? -halfomega - e : 0.0)

    for k=1:size(bands,2)
        col = real.(bands[:,k])
        sel = findall(inwindow, col)
        if length(sel) != dim #zone-edge degeneracy or truncation artifact: take the dim states closest to the window
            sel = partialsortperm(windowdist.(col), 1:dim)
        end
        Fbands[:,k] = sort!(fold.(col[sel]))
    end

    return Fbands
end

"""getfirstFBZ(bands::AbstractArray, M::Integer)
    Legacy method kept for backward compatibility: selects the middle dim rows of the energy-sorted
    bands. This is only correct as long as no quasienergy band crosses the edge of the Floquet
    Brillouin zone; otherwise a replica copy is included and the true band dropped. Pass the drive
    frequency omega (see the three-argument method) to get the correct windowed selection."""
function getfirstFBZ(bands::AbstractArray, M::Integer)
    dim = Int(size(bands,1)/(2*M+1)) #dimension of the Bloch Hamiltonian
    Fbands = zeros(dim,size(bands,2))

    for i=1:dim
        Fbands[i,:] = bands[M*dim+i,:]
    end

    return Fbands
end

"""transformeigvecs(U::AbstractMatrix, M::Integer, omega::Real, t::Number)
    Converts a matrix of eigenvectors U obtained from solving the infinite (truncated at M)
    Floquet problem with drive frequency omega to eigenvectors in time space (at time t)."""
function transformeigvecs(U::AbstractMatrix, M::Integer, omega::Real, t::Number)
    dim = Int(size(U,1)/(2*M+1)) #dimension of the Bloch Hamiltonian
    Unew = zeros(ComplexF64,dim, size(U,2))

    for i=1:size(U,2)
        Unew[:,i] = transform(U[:,i], dim, M, omega, t)
    end

    return Unew
end

"""Converts an eigenvector vec obtained from solving the infinite (truncated at M)
    Floquet problem with drive frequency omega to an eigenvector in time space (at time t)."""
function transform(vec::Array{ComplexF64,1}, dim::Integer, M::Integer, omega::Real, t::Number)
    vecnew = zeros(ComplexF64,dim)

    for i=0:(2*M)
        vecnew += exp(-im*(i-M)*omega*t).*vec[i*dim+1:(i+1)*dim]
    end

    return vecnew
end


import ..Structure.Paths
import ..Spectrum: dim

"""
    keepfirstFBZ!(data::Spectrum.BandData, H::FloquetOperator)

Throw away all band data except for the first Floquet Brillouin zone,
i.e., the quasienergies inside (-omega/2, omega/2] with omega the drive frequency.
"""
function keepfirstFBZ!(data::Spectrum.BandData, H::FloquetOperator)
    bands = getfirstFBZ(data.bands, H.M, H.drive.omega)
    data.bands = bands
    data
end

# Allow passing FloquetOperator directly to Spectrum.getbands by defining
# the appropriate dim methods here (no need for a wrapper function).
dim(HF::FloquetOperator, x::Number) = size(HF(x), 1)
dim(HF::FloquetOperator, x::AbstractVector) = size(HF(first(x)), 1)
dim(HF::FloquetOperator, x::AbstractMatrix) = size(HF(first(eachcol(x))), 1)


# """
#     geteigvecs(H, drive::periodicDrive, M::Integer; timespace=false)

# Returns a matrix of Floquet eigenstates of a system with time-averaged Hamiltonian H,
# driven by the periodic potential contained in drive.
# M denotes the frequency cutoff and therefore also the length of the eigenvectors in frequency space.
# For a d-dimensional Hamiltonian, the Fourier-space eigenvectors will have length d*(2*M+1).
# If timespace=true, these eigenvectors are converted back to time space and will have length dim.
# """
# function Spectrum.geteigvecs(HF::FloquetOperator; timespace=false)

#     HFloq = getFloquetHamiltonian(HF)
#     U = Spectrum.geteigvecs(HFloq)

#     if timespace
#         U = transformeigvecs(U, HF.M, 0)
#     end

#     return U
# end



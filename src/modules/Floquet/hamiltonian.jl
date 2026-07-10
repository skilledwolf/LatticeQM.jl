using LinearAlgebra: UniformScaling


"""
    getFloquetMatrix(M::Integer, H0::AbstractMatrix, drive::periodicDrive; hbar::Number=1.)

Build the effective Hamiltonian for several harmonic driving modes, i.e., \$V(t) = \\sum_j V_j \\exp(i*j*n_j*\\omega)\$.

M: truncation
H0: time averaged Hamiltonian (i.e., has already absorbed the constant m=0 mode), must be matrix
drive: periodicDrive object

The drive V(t) is assumed Hermitian, which implies \$V_{-n} = V_n^\\dagger\$ for its Fourier
components. Harmonics whose partner -n is not stored in the drive are completed
automatically with the adjoint, so the Floquet matrix is Hermitian even if only one member
of a ±n pair is given. If both +n and -n are stored they must be adjoints of each other
(checked within tolerance, otherwise an `ArgumentError` is thrown); the pair is then
symmetrized so the result is exactly Hermitian.
"""
function getFloquetMatrix(M::Integer, H0::AbstractMatrix, drive::periodicDrive; hbar::Number=1.)
    S = UniformScaling(drive.omega*hbar)
    d = size(H0,1) #dimension of the not-driven system
    H = spzeros(ComplexF64 ,d*(2*M+1),d*(2*M+1)) #store the effective Hamiltonian here try sparse arrays for now

    #collect the harmonics into one operator per n, summing duplicate entries of drive.ns
    harmonics = Dict{Int,Matrix{ComplexF64}}()
    for (j,n)=enumerate(drive.ns)
        n == 0 && continue #constant part belongs into H0 (see docstring)
        V = get!(() -> zeros(ComplexF64, d, d), harmonics, Int(n))
        V .+= drive.operators[j]
    end

    #Hermitian V(t) implies V_{-n} = V_n': complete missing partners, verify stored ones
    stored = sort!(collect(keys(harmonics)))
    tol = sqrt(eps(Float64))
    for n in stored
        n > 0 || continue
        if -n in stored
            if !isapprox(harmonics[-n], adjoint(harmonics[n]); rtol=tol, atol=tol)
                throw(ArgumentError("harmonics n=$(n) and n=$(-n) of the drive are not " *
                    "adjoints of each other; V(t) must be Hermitian, i.e. V_{-n} = V_n' " *
                    "(non-Hermitian drives are not supported)"))
            end
            Vsym = (harmonics[n] + adjoint(harmonics[-n]))/2 #symmetrize for exact Hermiticity
            harmonics[n] = Vsym
            harmonics[-n] = adjoint(Vsym)
        else
            harmonics[-n] = adjoint(harmonics[n])
        end
    end
    for n in stored #negative-only harmonics get their positive partner
        n < 0 && !(-n in stored) && (harmonics[-n] = adjoint(harmonics[n]))
    end

    for m=0:(2*M)
        i = m*d+1
        H[i:i+d-1, i:i+d-1] = H0 + (M-m)*S #fill the diagonal blocks

        for n in sort!(collect(keys(harmonics))) #iterate over the different drive frequencies
            if m < 2*M+1-abs(n) #off-diagonals have less blocks make sure to not run out of bounds
                if any(!iszero, harmonics[n]) #do not fill in matrices containing only zeros
                    if n>0
                        H[i:i+d-1, i+n*d:i+(n+1)*d-1] += harmonics[n] #fill the upper-right off-diagonal blocks
                    else
                        H[i-n*d:i+(-n+1)*d-1, i:i+d-1] += harmonics[n] #fill the lower-left off-diagonal blocks
                    end
                end
            end
        end
    end

    return dropzeros!(H) #get rid of 0 entries which werre introduced through H0 and operators
end

"""
    getFloquetHamiltonian(M::Number, H0, drive::periodicDrive)  
    getFloquetHamiltonian(M::Number, H0, drive::Function)

Build the effective Hamiltonian for several harmonic driving modes, i.e., \$V(t) = \\sum_j V_j \\exp(i*j*n_j*\\omega)\$.

M: truncation  
H0: time averaged Hamiltonian (i.e., has already absorbed the constant m=0 mode), must be callable
drive: periodicDrive object
"""

mutable struct FloquetOperator
    H0
    drive
    M::Integer
end

(HF::FloquetOperator)(k; kwargs...) = getFloquetMatrix(HF.M, HF.H0(k), HF.drive; kwargs...)

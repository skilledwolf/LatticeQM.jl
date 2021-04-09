using LinearAlgebra: UniformScaling


"""
    getFloquetHamiltonian(M::Integer, H0::AbstractMatrix, drive::periodicDrive; hbar::Number=1.)
    getFloquetHamiltonian(M::Number, H0::Function, drive::periodicDrive)  
    getFloquetHamiltonian(M::Number, H0::Function, drive::Function)

Build the effective Hamiltonian for several harmonic driving modes, i.e., \$V(t) = \\sum_j V_j \\exp(i*j*n_j*\\omega)\$.

M: truncation  
H0: time averaged Hamiltonian (i.e., has already absorbed the constant m=0 mode)  
drive: periodicDrive object
"""
function getFloquetHamiltonian(M::Integer, H0::AbstractMatrix, drive::periodicDrive; hbar::Number=1.)
    S = UniformScaling(drive.omega*hbar)
    d = size(H0,1) #dimension of the not-driven system
    H = spzeros(ComplexF64 ,d*(2*M+1),d*(2*M+1)) #store the effective Hamiltonian here try sparse arrays for now


    for m=0:(2*M)
        i = m*d+1
        H[i:i+d-1, i:i+d-1] = H0 + (M-m)*S #fill the diagonal blocks

        for (j,n)=enumerate(drive.ns) #iterate over the different drive frequencies
            if m < 2*M+1-abs(n) #off-diagonals have less blocks make sure to not run out of bounds
                if any(i->i!=0, drive.operators[j]) #do not fill in matrices containing only zeros
                    if n>0
                        H[i:i+d-1, i+n*d:i+(n+1)*d-1] += drive.operators[j] #fill the upper-right off-diagonal blocks
                    end
                    if n<0
                        H[i-n*d:i+(-n+1)*d-1, i:i+d-1] += drive.operators[j] #fill the lower.left off-diagonal blocks
                    end
                end
            end
        end
    end

    return dropzeros!(H) #get rid of 0 entries which werre introduced through H0 and operators
end

getFloquetHamiltonian(M::Number, H0::Function, drive::periodicDrive) = k -> getFloquetHamiltonian(M, H0(k), drive)
getFloquetHamiltonian(M::Number, H0::Function, drive::Function) = k -> getFloquetHamiltonian(M, H0(k), drive(k))


mutable struct FloquetOperator
    H0
    drive
    M::Integer
end

getFloquetHamiltonian(HF::FloquetOperator; kwargs...) = getFloquetHamiltonian(HF.M, HF.H0, HF.drive; kwargs...)
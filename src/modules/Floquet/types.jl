"""
    mutable struct periodicDrive

Type representing a time-periodic potential applied to a lattice.
Contains the fields:  
 - omega: Frequency of the potential.  
 - operators: An array of operators (Matrices) representing the Fourier components of the potential.  
 - ns: An array of integers containing the modes which are not zero. (Should not contain 0)  
    This implies that operators and ns must have the same length and must be ordered such that  
    operators[i] corresponds to the ns[i] fourier component for all i.

Example:  
    The following returns a drive corresponding to a periodic potentisl V(t):  
    V(t) = [0 1; -1 0]*exp(-im*2.5*t) + [0 -1; 1 0]*exp(3*im*2.5*t)


    drive = periodicDrive(2.5, [[0 1; -1 0], [0 -1; 1 0]], [-1,3])
"""
mutable struct periodicDrive
    omega::Float64 #circular frequency of the drive
    operators::AbstractArray{AbstractMatrix, 1} #List of operators corresponding to Fourier components of the drive
    ns::AbstractArray{Integer, 1} #list of multiples of omega to include in the drive
end

"""
    periodicDrive(omega::Number=1)

Return an empty periodic drive with frequency omega.
"""
function periodicDrive(omega::Number=1)
    periodicDrive(omega, [], [])
end

"""
    periodicDrive(omega::Number, funcFourier::AbstractArray{Number, 1}, operator::AbstractMatrix{Number})

Allows for initializing a periodicDrive using only one operator and Fourier components of a periodic function
denoting the amplitude (over time) with which the operator is applied.
"""
function periodicDrive(omega::Number, funcFourier::AbstractArray{Number, 1}, operator::AbstractMatrix{Number})
    M = convert(Integer, (length(funcFourier)-1)/2)
    ns = collect(-M:M)
    operators = [A.*operator for A in funcFourier]
    periodicDrive(omega, operators, ns)
end

"""
    addFreq!(drive::periodicDrive, operator::AbstractMatrix, n::Integer)

Add a single harmonic mode to an existing periodicDrive.
"""
function addFreq!(drive::periodicDrive, operator::AbstractMatrix, n::Integer)
    append!(drive.operators, [operator])
    append!(drive.ns, [n])
end


"""
    addcos!(drive::periodicDrive, operator::AbstractMatrix, n::Integer)

Add the operator applied with a cosine of frequency n*omega to an existing periodicDrive.
"""
function addcos!(drive::periodicDrive, operator::AbstractMatrix, n::Integer)
    append!(drive.operators, 0.5.*[adjoint(operator), operator])
    append!(drive.ns, [-n n])
end

function Base.show(io::IO, drive::periodicDrive)
    println(io, "Angular frequency of the periodic drive: ", drive.omega)
    println(io, "Non-zero harmonic modes: ", drive.ns)
    println(io, "Fourier components of the driving operator: ", drive.operators)
end

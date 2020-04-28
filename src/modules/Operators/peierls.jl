
"""
    peierls!(hops, lat, phase)

This method implements Peierl's substitution on a tight-binding model by
adding the appropriate phases to each hopping amplitude in the Hamiltonian
given by `hops`.

The phase function `phase` must be a function with the signature
    phase(r1::AbstractVector, r2::AbstractVector)

Note:
- this should be the last step when constructing a tight-binding 
  Hamiltonian
- You need to make sure that the phases that you add do not break
  translational symmetries. This is not a trivial matter and may
  lead to unexpected/undetected mistakes.
"""
function peierls!(hops, lat::Lattice, phase::Function)
    N = countorbitals(lat)
    D = hopdim(hops)
    d = div(D,N) # if spinhalf then d=2, if spinless d=1
    @assert d == 1 # only spinless (testing)
    @assert D == N*d #consistency check

    A = getA(lat)
    X = positions(lat)

    for (R, h) = hops
        δa = A * R
        for i=1:N,j=1:N
            hops[R][i,j] = hops[R][i,j] * exp(1im*2π*phase(X[:,i]+δa, X[:,j]))
            # hops[R][1+d*(i-1):d+d*(i-1),1+d*(j-1):d+d*(j-1)] .*= exp(1.0im*phase(X[:,i]+δa, X[:,j])) #sparse(I,d,d)
        end
    end

    hops
end


function peierls!(hops, lat::Lattice, B::AbstractVector)
    phase(args...) = uniformfieldphase(args...; B=B)
    peierls!(hops, lat, phase)
end


function peierlsinplane!(hops, lat::Lattice, B::AbstractVector)
    phase(args...) = uniformfieldphase_inplane(args...; B=B)
    peierls!(hops, lat, phase)
end


"""
    uniformfieldphase(r1,r2; B)

    The proper relative phase between lattice positions r1 and r2 in presence
    of magnetic field B=(B_1,B_2,...). Calculated in Coulomb gauge!

    The magnetic field should be in units of flux quanta \$\\phi_0=h/e\$.
"""
function uniformfieldphase(r1::T,r2::T; B::AbstractVector) where T<:AbstractVector
    @assert size(r1)==size(r2)==size(B) "Vectors must have same length"
    D = size(r1,1)
    @assert 1 < D < 4 "So far only 2D and 3D vectors are supported"

    cross2D(x, y) = [0, 0, x[1] * y[2] - x[2] * y[1]]
    f = (D==2) ? cross2D : cross

    dot(f((r1+r2)/2, (r2-r1)/2), B)
end


"""
    uniformfieldphase_inplane(r1,r2; B)

    The proper relative phase between lattice positions r1 and r2 in presence
    of magnetic field B=(B_1,B_2,0). 

    Here we use the gauge A = (B2 z, - B1 z, 0)

    The magnetic field should be in units of flux quanta \$\\phi_0=h/e\$.
"""
function uniformfieldphase_inplane(r1::T,r2::T; B::AbstractVector) where T<:AbstractVector
    @assert size(r1)==size(r2) "Vectors must have same length"
    @assert size(B,1)==2 "In-plane field must be 2D."
    D = size(r1,1)
    @assert 1 < D < 4 "So far only 2D and 3D vectors are supported"

    cross2D(x, y) = x[1] * y[2] - x[2] * y[1]

    if D<3
        z1=0;  z2=0
    else
        # z1=r1[3]+1.5;  z2=r2[3]+1.5
        z1=r1[3]+1.5;  z2=r2[3]+1.5
    end

    return cross2D(r2-r1, B) * (z1+z2)/2
end
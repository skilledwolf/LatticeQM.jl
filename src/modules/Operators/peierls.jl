import ..Structure.Lattices
import ..TightBinding

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
    N = Lattices.countorbitals(lat)
    D = TightBinding.hopdim(hops)
    d = div(D,N) # if spinhalf then d=2, if spinless d=1
    @assert d == 1 # only spinless (testing)
    @assert D == N*d #consistency check

    A = Lattices.getA(lat)
    X = Lattices.positions(lat)

    for (R, h) = hops
        δa = A * R
        for i=1:N,j=1:N
            hops[R][i,j] = hops[R][i,j] * exp(1im*2π*phase(X[:,i]+δa, X[:,j]))
            # hops[R][1+d*(i-1):d+d*(i-1),1+d*(j-1):d+d*(j-1)] .*= exp(1.0im*phase(X[:,i]+δa, X[:,j])) #sparse(I,d,d)
        end
    end

    hops
end

"""
    peierls!(hops, lat, B)

Add Peierls phases to operator `hops` on lattice geometry `lat` for the uniform
magnetic field `B=(B1,B2,B3)`. Uses Coulomb gauge, see uniformfieldphase(...).

"""
function peierls!(hops, lat::Lattice, B::AbstractVector)
    phase(args...) = uniformfieldphase(args...; B=B)
    peierls!(hops, lat, phase)
end

"""
    peierls!(hops, lat, B)

Add Peierls phases to operator `hops` on lattice geometry `lat` for the uniform
in-plane magnetic field `B=(B1,B2)`. Uses in-plane gauge, see uniformfieldphase_inplane(...).

"""
function peierlsinplane!(hops, lat::Lattice, B::AbstractVector)
    phase(args...) = uniformfieldphase_inplane(args...; B=B)
    peierls!(hops, lat, phase)
end


# function peierlsoutplane!(hops, lat::Lattice, phase::Function)
#     N = countorbitals(lat)
#     D = hopdim(hops)
#     d = div(D,N) # if spinhalf then d=2, if spinless d=1
#     @assert d == 1 # only spinless (testing)
#     @assert D == N*d #consistency check

#     # phasebar(r1,r2)= -(phase(r1,r2)-phase(r2,r1))/2
#     # phasebar(ri, rj, R) = phase(ri,rj)+(phase(ri+rj,R)-phase(R,ri+rj))/2

#     A = getA(lat)
#     X = positions(lat)

#     for (R, h) = hops
#         δa = A * R
#         for i=1:N,j=1:N
#             # hops[R][i,j] *= exp(1im*2π*phasebar(X[:,i]+δa, X[:,i]+X[:,j])) #hops[R][i,j] * exp(1im*2π*phasebar(X[:,i]+δa, X[:,i]+X[:,j]))
#             hops[R][i,j] *= exp(1im*2π*phase(X[:,i], X[:,j], δa)) #hops[R][i,j] * exp(1im*2π*phasebar(X[:,i]+δa, X[:,i]+X[:,j]))
#         end
#     end

#     hops
# end

function peierlsoutplane!(hops, lat::Lattice, phase::Function)
    N = Lattices.countorbitals(lat)
    D = TightBinding.hopdim(hops)
    d = div(D,N) # if spinhalf then d=2, if spinless d=1
    @assert d == 1 "Only spin-less implemented so far for peierl phase." # only spinless (for now)

    A = Lattices.getA(lat)
    X = Lattices.positions(lat)

    ϕ(i,j,R) = phase(X[:,i], X[:,j], A*R)

    peierlsoutplane!(hops, ϕ)
end

function peierlsoutplane!(hops, phase::Function)
    D = TightBinding.hopdim(hops)

    for (R, h) = hops
        for i=1:D, j=1:D
           hops[R][i,j] *= exp(1im*2π*phase(i,j,R)) #hops[R][i,j] * exp(1im*2π*phasebar(X[:,i]+δa, X[:,i]+X[:,j]))
        end
    end

    hops
end


import ..Structure.Lattices

"""
    peierlsoutplane(hops, lat, p, q)

Add uniform out-of-plane magnetic field B=(0,0,B3) to operator `hops` on lattice geometry `lat`,
such that the flux per unit cell is Φ=p/q.

This method automatically constructs and returns the correct superoperartor and magnetic supercell.
"""
function peierlsoutplane(hops, lat::Lattice, p::Int, q::Int; verbose=true)
    @assert gcd(p,q)<2 "Integers p,q for flux Φ=p/q should be coprime. (here p=$p, q=$q)"
    verbose ? print("Building magnetic supercell for Φ=p/q=$p/$q.") : nothing

    Θ = uniformfieldphase_outplane(lat, p/q)
    Φ = (r,R) -> Θ((r+R)/2,r-R)

    A = Lattices.getA(lat)

    mhops = deepcopy(hops)
    peierlsoutplane!(mhops, lat, (ri,rj,R)->Φ(ri,rj)+(Θ(ri+rj,R)-Θ(R,ri+rj))/2) # add the cell-position independent phases first

    mlat  = Lattices.superlattice(lat, [1,q]) # make correct supercell
    # mhops = TightBinding.superlattice(mhops, [1,q], (r,R)->exp(1im*2π*(r[2]*R[1])*(p/q)) ) # make correct supercell hamiltonian
    mhops = superlattice(mhops, [1,q], (r,R)-> exp(1im*2π*Φ(A*r,A*R)) )

    mhops, mlat
end


using ..Spectrum: bandmatrix, getdos!
using ProgressMeter

"""
    hofstadter(hops, lat, Q)

Determines the energie spectrum as function of rational magnetic flux \$\\Phi=p/q\$,
where p, q are coprime integers with 1<= q <= Q and 1<=p<q.
Returns a list of fluxes and a list of energies at each flux.

In this implementation we evaluate at the \$\\Gamma\$-point of the magnetic cell.
"""
function hofstadter(hops, lat::Lattice, Q::Int; kwargs...)

    fluxes = [(p,q) for q=1:Q for p=1:q-1 if gcd(p,q)<2]
    energies = Vector{Float64}[]

    k0 = zeros(Float64, Lattices.latticedim(lat), 1)

    @showprogress 2 "Iterating through flux... " for (p,q)=fluxes
        mhops, mlat = Operators.peierlsoutplane(hops, lat, p, q; verbose=false)
        append!(energies, [vec(bandmatrix(mhops, k0; kwargs...)[1])])
    end

    fluxes = map(x->x[1]//x[2], fluxes)
    p = sortperm(fluxes)
    fluxes = fluxes[p]
    energies = energies[p]

    fluxes, energies
end

"""
    hofstadter_dos(hops, lat, q_max::Int, (f_min, f_max), N=300)

See hofstadter_dos(hops, lat, q_max, frequencies).

Returns a list of fluxes, frequencies and a the dos at each flux.
"""
function hofstadter_dos(hops, lat::Lattice, Q::Int, ω::Tuple, N::Int=300; kwargs...)
    ωs = LinRange(ω...,N)
    fluxes, DOS = hofstadter_dos(hops, lat, Q, ωs; kwargs...)
    fluxes, ωs, DOS
end

"""
    hofstadter_dos(hops, lat, q_max::Int, frequencies::AbstractVector; klin=100, Γ=0.05)

Determines the energie spectrum as function of rational magnetic flux \$\\Phi=p/q\$,
where p, q are coprime integers with 1<= q <= q_max and 1<=p<q.
Note that q_max determines the size of the largest magnetic supercell.

The density of states (DOS) at each flux is calculated on a discrete k-grid (resolution give by `klin`).
The parameter Γ determines the energy broadening when calculating DOS.

Returns a list of fluxes and a the dos at given frequencies at each flux.
"""
function hofstadter_dos(hops, lat::Lattice, Q::Int, frequencies::AbstractVector; klin::Int=100, Γ::Real=0.05)

    fluxes = [(p,q) for q=1:Q for p=1:q-1 if gcd(p,q)<2]
    DOS = zeros(length(frequencies), length(fluxes))

    @showprogress 2 "Iterating through flux... " for (i,(p,q))=enumerate(fluxes)
        mhops, mlat = Operators.peierlsoutplane(hops, lat, p,q; verbose=false)
        
        getdos!(view(DOS, :, i), mhops, frequencies; klin=round(Int,klin/√q), Γ=Γ)
    end

    fluxes = map(x->x[1]//x[2], fluxes)
    p = sortperm(fluxes)

    fluxes[p], DOS[:,p]
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

    dot(f(r1, r2), B/2)
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
    # end
    else
        ### z1=r1[3]+1.5;  z2=r2[3]+1.5
        z1=r1[3];  z2=r2[3]
    end
    # z1=r1[3]; z2=r2[3]

    return cross2D(r2-r1, B) * (z1+z2)/2
end

import ..Structure.Lattices

"""
    uniformfieldphase_outplane(lat, Φ)

Passes the lattice vectors a1 and a2 to uniformfieldphase_outplane(a1,a2,Φ).
"""
function uniformfieldphase_outplane(lat::Lattice, args...; kwargs...)
    @assert Lattices.latticedim(lat)==2 "Only implemented for 2D lattices."
    @assert Lattices.spacedim(lat)>2 "Only implemented for lattices in at least 3D space."
    A = Lattices.getA(lat)
    uniformfieldphase_outplane(A[:,1], A[:,2], args...; kwargs...)
end

"""
    uniformfieldphase_outplane(a1,a2)

This returns the phase function that respects translational symmetry along \$a_1\$
but not along \$a_2\$. For rational flux Φ=p/q, this gauge is periodic in a2'=q*a2.

The magnetic field should be in units of flux quanta \$\\phi_0=h/e\$.
"""
function uniformfieldphase_outplane(a1::T,a2::T, Φ::Number) where T<:AbstractVector
    @assert size(a1)==size(a2) "Vectors a1, a2 must have same length"
    @assert size(a1,1) == 3 "So far only 3D vectors are supported"

    cross2D(x, y) = x[1] * y[2] - x[2] * y[1]

    na1=norm(a1); na2=norm(a2)
    A = norm(cross(a1[1:3],a2[1:3]))
    e1 = a1/na1; e2 = a2/na2 # normalize

    zv = cross(e1[1:3], e2[1:3])
    s = norm(zv)
    zv = zv/s

    e1p = deepcopy(a1); e1p[1:3] .= cross(zv,e1[1:3])
    e2p = deepcopy(a2); e2p[1:3] .= cross(zv,e2[1:3])
    
    """
        phase(rM, δr)

        rM = (r1+r2)/2
        δr = r1-r2
    """
    function phase(rM::T,δr::T)
        Φ * dot(rM, e1p) * dot(δr,e2p) / (A*s)
    end

    return phase
end




# function peierlsoutplane2(hops, lat::Lattice, p::Int, q::Int; verbose=true)
#     @assert gcd(p,q)<2 "Integers p,q for flux Φ=p/q should be coprime. (here p=$p, q=$q)"
#     verbose ? print("Building magnetic supercell for Φ=p/q=$p/$q.") : nothing

#     mlat  = Structure.Lattices.superlattice(lat, [1,q])
#     mhops = TightBinding.superlattice(hops, [1,q])

#     phasefunc(r1,r2) = uniformfieldphase(r1, r2; B=[0,0,p/q])

#     peierls!(mhops, mlat, phasefunc)

#     mhops, mlat
# end

# function hofstadter2(hops, lat::Lattice, Q::Int)

#     lat0  = Structure.Lattices.superlattice(lat, [1,1])
#     hops0 = TightBinding.superlattice(hops, [1,-1])

#     fluxes = [(p,q) for q=1:Q for p=1:q-1 if gcd(p,q)<2]
#     energies = Vector{Float64}[]

#     k0 = zeros(Float64, Structure.latticedim(lat), 1)

#     @showprogress 2 "Iterating through flux... " for (p,q)=fluxes
#         mhops, mlat = Operators.peierlsoutplane2(hops0, lat0, p,q; verbose=false)
#         append!(energies, [vec(bandmatrix(mhops, k0)[1])])
#     end

#     fluxes = map(x->x[1]//x[2], fluxes)
#     p = sortperm(fluxes)
#     fluxes = fluxes[p]
#     energies = energies[p]

#     fluxes, energies
# end
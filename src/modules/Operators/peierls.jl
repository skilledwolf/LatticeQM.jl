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


using ..Spectrum: bandmatrix, getdos!, gapcherns, spectralgaps
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
    hofstadter_cherns(hops, lat, Q::Int, NX::Int=18; gaptol=1e-6, kwargs...)

Topological Hofstadter (Wannier) diagram: the Chern number of every spectral
gap as a function of rational magnetic flux \$\\Phi=p/q\$ per unit cell (coprime
p, q with 1 ≤ q ≤ Q and 1 ≤ p < q).

For each flux the magnetic supercell Hamiltonian is built with
[`peierlsoutplane`](@ref) and its gaps are labelled by [`Spectrum.gapcherns`](@ref) —
the cumulative non-abelian Chern number of the magnetic subbands below each gap
(= Hall conductance σxy in e²/h). The magnetic BZ is sampled on an `NX × NX`
grid; `NX ≳ 18` is recommended, as the densely packed subbands need a fine grid
to converge.

Filling is reported per *original* unit cell, `ν = n/q` for `n` subbands
filled, so the gaps fall on Diophantine/Středa lines `ν = s − C·Φ` with integer
Chern number `C` and band offset `s` — the sign is fixed by the flux sign of
[`peierlsoutplane`](@ref) together with the Berry orientation of
[`Spectrum.gapcherns`](@ref) and is pinned by a regression test
(honeycomb, all gaps at q ≤ 5).

Returns a `NamedTuple` of equal-length vectors `(; flux, filling, energy, width, chern)`,
one entry per (flux, gap): `flux::Rational` is Φ=p//q, `filling` is ν per unit
cell, `energy` is the gap-centre energy, `width` is the gap size, and `chern`
is the gap Chern number. With the default near-zero `gaptol` the list includes
numerically tiny pseudo-gaps; filter or weight by `width` for a clean diagram.

`kwargs` are forwarded to the diagonalisation (`multimode`, `executor`, `format`).
"""
function hofstadter_cherns(hops, lat::Lattice, Q::Int, NX::Int=18; gaptol::Real=1e-6, kwargs...)
    fluxlist = [(p, q) for q = 1:Q for p = 1:q-1 if gcd(p, q) < 2]

    flux    = Rational{Int}[]
    filling = Float64[]
    energy  = Float64[]
    width   = Float64[]
    chern   = Int[]

    @showprogress 2 "Hofstadter Chern numbers... " for (p, q) in fluxlist
        mhops, _ = peierlsoutplane(hops, lat, p, q; verbose=false)
        for g in gapcherns(mhops, NX; gaptol=gaptol, kwargs...)
            push!(flux,    p // q)
            push!(filling, g.n / q)
            push!(energy,  (g.elo + g.ehi) / 2)
            push!(width,   g.ehi - g.elo)
            push!(chern,   g.chern)
        end
    end

    (; flux, filling, energy, width, chern)
end


"""
    diophantine_cherns(p, q, r; Nw, smargin=0) -> Vector{Tuple{Int,Int}}

Candidate `(s, C)` labels of a Hofstadter gap from the TKNN / Středa–Wannier
Diophantine equation

    r = q·s − p·C,

for flux Φ=p/q (coprime) with `r` magnetic subbands filled below the gap and
`Nw` orbitals per primitive cell. `C ∈ ℤ` is the gap Chern number as reported
by [`Spectrum.gapcherns`](@ref) / [`hofstadter_cherns`](@ref) (the sign
convention those routines produce), and `s ∈ ℤ` the band-filling label; the
integrated density per primitive cell is `n = r/q = s − C·Φ` (in the
conventional Wannier-diagram form `n = s + t·Φ` the Středa slope is `t = −C`).

The equation fixes `C ≡ −p⁻¹ r (mod q)` only, so every branch inside the
heuristic window `−smargin ≤ s ≤ Nw + smargin` is returned, ordered by `|C|`.

!!! warning "The window and the |C| ordering are heuristics"
    The window bounds the ν-axis *intercept* of the gap's Středa line, not its
    filling: a gap family that exists only in a finite flux window can have a
    true intercept outside `[0, Nw]`, in which case the true branch is **not**
    among the `smargin = 0` candidates. This already happens for plain
    NN honeycomb: the gap at Φ = 4/5, r = 3 has (s, C) = (3, 3), outside
    `0 ≤ s ≤ Nw = 2` — `smargin ≥ 1` recovers it (pinned by a test).
    Likewise the smallest-`|C|` candidate is systematically wrong for the
    higher lines of a Dirac Landau fan (honeycomb near half filling has true
    labels C = ∓1, ∓3, ∓5, … at Φ → 0⁺, and once `|C| > q/2` an in-window
    alternative with smaller `|C|` exists in the same Diophantine class).
    A single candidate is an unambiguous label only up to these caveats;
    several means the branch must be pinned otherwise (a gap-family fit, or an
    exact [`Spectrum.gapcherns`](@ref) anchor) — do **not** blindly take the
    smallest `|C|` for a multi-orbital or Dirac-type model.
"""
function diophantine_cherns(p::Int, q::Int, r::Int; Nw::Int, smargin::Int=0)
    q == 1 && return [(r, 0)]                 # integer flux: C undefined, n = r
    invp = invmod(mod(p, q), q)
    C0 = mod(-invp * r, q)                    # C ≡ −p⁻¹ r (mod q), in 0:q-1
    out = Tuple{Int,Int}[]
    # s = (r + p·C)/q increases by p as C increases by q; sweep enough branches
    # to bracket the window on either side.
    lmax = cld(abs(r) + abs(p) * q + Nw + 2 * smargin + q, q) + 2
    for l in -lmax:lmax
        C = C0 + l * q
        s, rem = divrem(r + p * C, q)
        (rem == 0 && -smargin <= s <= Nw + smargin) && push!(out, (s, C))
    end
    sort!(out; by = sc -> (abs(sc[2]), sc[2]))
    out
end

"""
    hofstadter_gaplabels(hops, lat, Q::Int, NX::Int=18; gaptol=0.05, Nw=hopdim(hops), smargin=0)

Cheap (energies-only) Hofstadter gap labelling: the Diophantine/Středa Chern
label of every spectral gap as a function of rational flux Φ=p/q per unit cell
(coprime p, q with 1 ≤ q ≤ Q and 1 ≤ p < q).

For each flux the magnetic supercell is built with [`peierlsoutplane`](@ref) and
its gaps located by [`Spectrum.spectralgaps`](@ref) (no eigenvectors, no Berry —
much cheaper than [`hofstadter_cherns`](@ref)). Each gap is then labelled with
[`diophantine_cherns`](@ref): `chern` is the smallest-`|C|` candidate and `nsol`
the number of branches inside the `−smargin ≤ s ≤ Nw + smargin` window. `Nw` is
the number of orbitals per primitive cell. Every gap is reported: `nsol == 0`
means no in-window candidate exists (then `chern` is set to 0 and is
meaningless — widen `smargin`), `nsol > 1` flags a gap whose label needs
disambiguation (a gap-family fit or a [`hofstadter_cherns`](@ref) anchor).

The labels are Level-1 *guesses*; see the [`diophantine_cherns`](@ref) warning
for the systematic failure modes (out-of-window intercepts, Dirac Landau fans
— both occur already for plain NN honeycomb). Treat `chern` as a candidate to
be pinned by anchors, not as a result.

Filling is reported per primitive cell, `ν = r/q`. Returns a `NamedTuple` of
equal-length vectors `(; flux, filling, energy, width, q, r, chern, nsol)`.

This is the cheap Level-1 layer of the standard Hofstadter workflow; use
[`hofstadter_cherns`](@ref) for the exact Berry-curvature anchors.
"""
function hofstadter_gaplabels(hops, lat::Lattice, Q::Int, NX::Int=18;
                              gaptol::Real=0.05, Nw::Int=TightBinding.hopdim(hops),
                              smargin::Int=0)
    fluxlist = [(p, q) for q = 1:Q for p = 1:q-1 if gcd(p, q) < 2]

    flux    = Rational{Int}[]
    filling = Float64[]
    energy  = Float64[]
    width   = Float64[]
    qs      = Int[]
    rs      = Int[]
    chern   = Int[]
    nsol    = Int[]

    @showprogress 2 "Hofstadter gap labels... " for (p, q) in fluxlist
        mhops, _ = peierlsoutplane(hops, lat, p, q; verbose=false)
        for g in spectralgaps(mhops, NX; gaptol=gaptol)
            cands = diophantine_cherns(p, q, g.n; Nw=Nw, smargin=smargin)
            C = isempty(cands) ? 0 : first(cands)[2]  # smallest |C| (Level-1 guess)
            push!(flux,    p // q)
            push!(filling, g.n / q)
            push!(energy,  (g.elo + g.ehi) / 2)
            push!(width,   g.ehi - g.elo)
            push!(qs,      q)
            push!(rs,      g.n)
            push!(chern,   C)
            push!(nsol,    length(cands))
        end
    end

    (; flux, filling, energy, width, q=qs, r=rs, chern, nsol)
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
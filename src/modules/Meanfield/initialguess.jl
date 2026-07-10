using ..TightBinding
using ..Operators: getoperator
using ..TightBinding: zerokey
import ..Utils

"""
    spindensitymatrix(d::Vector=[1,0,0]) --> Matrix{Float64}
    spindensitymatrix(x::Real=1) --> Matrix{Float64}
    spindensitymatrix(s::Symbol) --> Matrix{Float64}

Returns a spin-1/2 density matrix for given spin orientiation d=[dx,dy,dz].
Instead, one can also specify the configuration symbolically, in which
case s should be one of `:up,:down,:upx,:downx,:upy,:downy`.

If vector has norm larger one it gets normalized, otherwise it is left unscaled.

"""
spindensitymatrix(x::Real=1) = spindensitymatrix([0,0,x])
function spindensitymatrix(args...; kwargs...)
    M = zeros(ComplexF64, 2,2)
    spindensitymatrix!(M, args...; kwargs...)
end

function spindensitymatrix!(M::AbstractMatrix, d::Vector=[1,0,0]; i=1)
    @assert length(d)<=3 "Spin direction must be length 3."
    d0 = norm(d)
    if d0 > 1
        d = d/d0
    end

    δ = 2*(i-1)
    M[1+δ:2+δ,1+δ:2+δ] .= 0.5 .* (Utils.σ0 .+ sum(d[i_] .* Utils.σs[i_] for i_=1:length(d)))

    M
end

function spindensitymatrix!(M::AbstractMatrix, s::Symbol; kwargs...)
    if s==:left || s==:upx || s==:x
        d = [1,0,0]
    elseif s==:right || s==:downx
        d = [-1,0,0]
    elseif s==:front || s==:upy || s==:y
        d = [0,1,0]
    elseif s==:back || s==:downy
        d = [0,-1,0]
    elseif s==:up || s==:z
        d = [0,0,1]
    elseif s==:down
        d = [0,0,-1]
    elseif s==:random || s==:rand
        θ = π*rand()
        ϕ = 2π*rand()
        d = [sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)]
    elseif s==:randomXY || s==:randXY
        ϕ = 2π * rand()
        d = [cos(ϕ), sin(ϕ), 0]
    elseif s==:randomZ || s==:randZ
        d = [0,0,(-1)^rand(Bool)]
    end

    spindensitymatrix!(M, d; kwargs...)
end


"""
    mapspindensitymatrix(vs::AbstractVector, Is::AbstractVector{Int}, N::Int) --> SparseMatrix{Complex}
    mapspindensitymatrix(vs::AbstractVector, Is::AbstractVector{Int}) = mapspindensitymatrix(vs, Is, length(Is)) --> SparseMatrix{Complex}
    mapspindensitymatrix(vs::AbstractVector) --> SparseMatrix{Complex}

From a list of spin orientations vs (and optionally lattice indices Is, and optionally total number of lattice sites N),
a spin density matrix for multiple lattice sites is generated.
"""
mapspindensitymatrix(vs::AbstractVector) = mapspindensitymatrix(vs, 1:length(vs))
mapspindensitymatrix(v, N::Int) = mapspindensitymatrix(collect(Iterators.repeated(v, N)))
mapspindensitymatrix(vs::AbstractVector, Is::AbstractVector{Int}) = mapspindensitymatrix(vs, Is, length(Is))
function mapspindensitymatrix(vs::AbstractVector, Is::AbstractVector{Int}, N::Int)
    @assert length(Is)==length(vs) "List of indices Is must have same length as list of orientations vs."
    @assert N >= length(Is) "Target size must be at least as large as specified indices."

    M = spzeros(ComplexF64, 2*N, 2*N) # empty matrix
    for (i_,v) in zip(Is, vs) # fill up the matrix along block diagonals
        spindensitymatrix!(M, v; i=i_)
    end

    M
end



using ..Structure
import ..Structure.Lattices

function spinspiralangles(lat::Lattice, superperiods)

    R = eachcol(Lattices.positions(lat))
    slat = Lattices.superlattice(lat, superperiods)
    B = Lattices.getB(slat)

    2π * vec(sum(transpose(B)*R; dims=1))
end

import ..Structure

"""
    spinspiraldensitymatrix(lat::Lattice, superperiods; n::Vector=[0,0,1], v0::Vector=[1,0,0])

Create a spin spiral density matrix with periodicity given by a superlattice characterized by `superperiods' (see Structure.superlattice for details).
The rotation goes around the fixed axis `n=[n1,n2,n3]`. The vector `v0` specifies the spin direction in one site.

Useful to create initial states for mean field calculations.
"""
function spinspiraldensitymatrix(lat::Lattice, superperiods; n::Vector=[0,0,1], v0::Vector=[1,0,0])
    @assert length(n) == 3 "Normal vector for rotation must have length 3."
    @assert length(v0) == 3 "Initial vector for rotation must have length 3."

    v(θ) = Structure.rotation3D(θ,n)*v0
    θs = spinspiralangles(lat, superperiods)

    mapspindensitymatrix(v(θ) for θ=θs)
end


function setrandom!(ρ::Hops, kind=:nonlocal; filling=nothing)
    N = hopdim(ρ); @assert mod(N,2)==0 "The hopping matrix must have even dimension (i.e. spinful)."
    n = div(N,2)

    if kind==:local
        M = mapspindensitymatrix(:random, n)
        setzero!(ρ, M)
    elseif kind==:nonlocal
        for L in keys(ρ)
            ρ[L] .= 0
        end
        for L in keys(ρ)
            M = rand(ComplexF64, N, N); #M = (M+M')/2
            ρ[L] .+= M
            ρ[-L] .+= M'
        end
    elseif kind==:XY || kind==:xy
        M = mapspindensitymatrix(:randomXY, n)
        setzero!(ρ, M)
    elseif kind==:Z || kind==:z
        M = mapspindensitymatrix(:randomZ, n)
        setzero!(ρ, M)
    else
        error("Unrecognized request for mode '$kind'.")
    end

    # Project the diagonal of the zero-key block so tr(ρ[0]) = filling*N.
    # Without this, e.g. `:nonlocal` leaves a uniformly-distributed random
    # density (mean diagonal ≈ 1) that bears no relation to the target
    # filling — SCF then has to drift μ from a meaningless starting point.
    if filling !== nothing
        zk = zerokey(ρ)
        H0 = ρ[zk]
        target = filling * N
        actual = real(tr(H0))
        δ = (target - actual) / N
        @inbounds for i in 1:N
            H0[i, i] += δ
        end
    end

    ρ
end

function setferro!(ρ::Hops, d=:up) 
    N = hopdim(ρ); @assert mod(N,2)==0 "The hopping matrix must have even dimension (i.e. spinful)."
    n = div(N,2)

    if d==:random # makes sure that random only computes a single random vector!
        θ = π*rand()
        ϕ = 2π*rand()
        d = [sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)]
    elseif d==:randomXY
        ϕ = 2π * rand()
        d = [cos(ϕ), sin(ϕ), 0]
    end

    setzero!(ρ, mapspindensitymatrix(d, n))
end

function setantiferro!(ρ::Hops, lat::Lattice, d=:up)
    N = hopdim(ρ); @assert mod(N,2)==0 "The hopping matrix must have even dimension (i.e. spinful)."
    n = div(N,2)

    sublA, sublB = getoperator(lat, ["sublatticeA", "sublatticeB"])
    ρup = spindensitymatrix(d)

    setzero!(ρ, kron(sublA, ρup) + kron(sublB, Utils.σ0-ρup))
end

function initialguess(v::Hops, mode=:random, args...; lat=:nothing, kwargs...)
    N = hopdim(v); @assert mod(N,2)==0 "The hopping matrix must have even dimension (i.e. spinful)."
    n = div(N,2)

    # Density matrices are dense and complex regardless of v's backing
    # (v is often real/sparse). Build the accumulator explicitly instead of
    # relying on `zero(v)` — the dense-complex `Base.zero(::Hops)` overload
    # that used to provide this pirated the type-preserving TightBinding one.
    ρs = Hops(Dict(δL => zeros(ComplexF64, size(v[δL])) for δL in keys(v)))

    if mode==:random
        setrandom!(ρs, args...; kwargs...)
    elseif mode==:ferro
        setferro!(ρs, args...; kwargs...)
    elseif mode==:antiferro
        setantiferro!(ρs, lat, args...; kwargs...)
    elseif mode==:spiral
        setzero!(ρs, spinspiraldensitymatrix(lat, args...; kwargs...))
    else
        error("Unrecognized mode '$mode' in initialguess(...).")
    end

    ρs
end

"""
    enrichkeys!(ρ::Hops, klin::Int) → Hops

Expand `keys(ρ)` to cover the full Wigner-Seitz cell of the inverse k-grid
that `regulargrid(nk=klin^d)` produces, so the SCF can populate every
Fourier mode of `ρ_full(k)`. Keys are added with zero matrix blocks; any
existing block content is preserved.

This is the prerequisite for using `residual_norm=:commutator` as a hard
convergence criterion: with sparse `keys(ρ)` (the default for e.g. on-site
Hubbard `v`), `bloch!(ρ, k)` reconstructs only a low-rank approximation of
`ρ_full(k)`, and `‖[H, ρ]‖_F` hits a non-vanishing truncation floor.

# klin parity

`klin` should be **odd**. The WS cell `{-Lmax, …, Lmax}^d` for odd `klin`
is closed under `L → -L`, so the lattice-Hops Hermiticity condition
`ρ[L]^† = ρ[-L]` is satisfied without adding redundant keys. With even
`klin`, `hermitianize!` adds keys at the BZ boundary that double-count the
`k = ±½` Fourier mode, leaving a small (~10⁻²) residual floor.

# Cost

Memory: `klin^d × N²` complex doubles. For `klin = 9`, `d = 2`, `N = 4`
(graphene-Hubbard): 81 × 16 ≈ 1300 doubles. Negligible. For `klin = 9`,
`N = 1588` (TBG-N=11): ~ 200 MB — non-trivial but tractable.

# Example

```julia
ρ_init = Meanfield.initialguess(v, :antiferro; lat=lat)
Meanfield.enrichkeys!(ρ_init, 9)  # commutator residual now valid
ρ, ϵ, mf, conv, res = Meanfield.solvehartreefock(
    hops, v, ρ_init, 0.5; klin=9, residual_norm=:commutator, tol=1e-7)
# converges in ~8 iterations to res ≈ 1e-9
```
"""
function enrichkeys!(ρ::Hops, klin::Int)
    isodd(klin) || @warn "enrichkeys!: klin=$klin is even. The Hops " *
                         "Hermiticity check will add duplicate keys at the " *
                         "BZ boundary, leaving a small commutator residual " *
                         "floor. Use odd klin for machine-precision convergence."
    N = hopdim(ρ)
    isempty(ρ.data) && error("enrichkeys!: ρ must have at least one key " *
                              "to infer matrix type and lattice dimension")
    sample_key = first(keys(ρ.data))
    sample_mat = ρ.data[sample_key]
    d = length(sample_key)
    Lmax = klin ÷ 2
    for L_indices in Iterators.product(ntuple(_ -> -Lmax:Lmax, d)...)
        L = typeof(sample_key)(collect(L_indices))
        if !haskey(ρ.data, L)
            ρ.data[L] = zero(sample_mat)
        end
    end
    ρ
end
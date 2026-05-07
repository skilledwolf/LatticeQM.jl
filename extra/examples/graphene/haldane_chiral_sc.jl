using LinearAlgebra
using Printf
using Random
using Plots

using LatticeQM

# ---------------------------------------------------------------------------
# Haldane model with on-site sublattice mass Δ, onsite repulsion U, and
# attractive NN interaction V.
#
# Step 1: build a spinful Haldane Hamiltonian with on-site sublattice mass
#         (−Δ on A, +Δ on B); plot the valley-resolved bands. With Δ chosen
#         so the K- and K'-valley masses |Δ ± m_H| are very different, one
#         valley sits at low energy and the other is pushed up — the
#         "valley-polarised" low-energy band structure.
# Step 2: run a BdG mean-field solve with onsite U > 0 and V_NN < 0. The
#         onsite repulsion suppresses the extended-s basin; with the full
#         density Hartree-Fock channel retained, a chiral d±id NN singlet seed
#         converges to a finite m = +1 solution for the parameters below.
# Step 3: read off the singlet pairing amplitude on each of the three NN
#         A→B bonds, plot them in the complex plane, and check the C3
#         winding number.
#
# Geometry. The three A→B nearest-neighbour bonds are indexed by the unit-
# cell offset δR carrying B (verified empirically against
# `Lattices.getneighbors(lat, 1.0)`):
#   δR = ( 0, 0):  δ₁ = (+1,  0)
#   δR = (−1, 0):  δ₂ = (−1/2, +√3/2)
#   δR = ( 0,−1):  δ₃ = (−1/2, −√3/2)
# The C3 rotation R(2π/3) cyclically permutes them: δ₁ → δ₂ → δ₃ → δ₁.
# A pairing channel of angular momentum m therefore satisfies
#   arg Δ(δ_{a+1}) − arg Δ(δ_a) = 2π m / 3   (mod 2π).
# m = 0   → extended s-wave (or sublattice-anisotropic s)
# m = ±1  → chiral d±id (or p±ip, depending on the spin channel)
# ---------------------------------------------------------------------------

Random.seed!(2026)

# Model parameters
const t2_val   = 0.10    # Haldane NNN amplitude (Haldane mass at K = 3√3 t₂ ≈ 0.520)
const ϕ_val    = π/2     # Haldane phase (opens the topological gap)
const Δ_mass   = 0.40    # On-site sublattice mass: ±Δ_mass on A/B
                         #   K-valley gap  = 2|Δ − m_H| ≈ 0.24  ← low-energy valley
                         #   K'-valley gap = 2|Δ + m_H| ≈ 1.84
const U_onsite = 1.0     # Onsite Hubbard repulsion U n↑n↓
const V_attr   = -2.5    # Attractive nearest-neighbour interaction
const filling  = 0.55    # Slight doping just above charge neutrality. With
                         # the K-valley upper band starting at +0.12 above
                         # the Dirac point, μ ≈ +0.34 here puts a small Fermi
                         # pocket entirely in the K band — single-valley
                         # carriers that pair via the NN attraction. Heavier
                         # doping (filling ≳ 0.6) drives the SCF to a
                         # C3-broken nematic; charge neutrality (filling ≡
                         # 0.5) gives BCS too small for finite-tol SCF to
                         # distinguish from Δ ≡ 0.
const klin_scf = 25      # Odd → BZ closed under k → −k (clean Hops Hermiticity)
const T_scf    = 0.01

# ---------------------------------------------------------------------------
# Step 1 — build the single-particle Hamiltonian and plot valley-resolved bands
# ---------------------------------------------------------------------------
println(":: Step 1 — single-particle Haldane + Δ Hamiltonian")

lat     = Geometries.honeycomb()
ks_path = kpath(lat; num_points=200)

# `addhaldane!` writes only into the (1,1)-block per orbital (it predates spin
# support); build it on a spinless Hops first, then promote with `addspin`
# so the Haldane phase appears identically on both spin sectors.
h_sp = Operators.graphene(lat; mode=:nospin)
Operators.addhaldane!(h_sp, lat, t2_val; ϕ=ϕ_val)
# `addsublatticeimbalance!(h, lat, x)` puts (x/2)·(sublattice − 1/2) on the
# diagonal, i.e. −x/2 on A and +x/2 on B; pass 2·Δ_mass to land at ±Δ_mass.
Operators.addsublatticeimbalance!(h_sp, lat, 2 * Δ_mass)
h0 = TightBinding.addspin(h_sp, :spinhalf)

valley_op = Operators.valley(lat; spinhalf=true)
bands     = getbands(h0, ks_path, valley_op)

p_bands = plot(bands;
               ylabel = "\$\\varepsilon/t\$",
               colorbar_title = "valley",
               markersize     = 2.0,
               markercolor    = :PiYG,
               size           = (380, 300),
               title          = @sprintf("Haldane: Δ=%.2f, t₂=%.2f, ϕ=π/2",
                                          Δ_mass, t2_val))

# ---------------------------------------------------------------------------
# Step 2 — BdG mean-field SCF with onsite U and attractive NN interaction
# ---------------------------------------------------------------------------
println(":: Step 2 — BdG SCF, U = $(U_onsite), V_NN = $(V_attr), filling = $(filling)")

hbdg = BdGOperator(h0)
v    = Operators.getshortrangedpotential(lat, U_onsite, V_attr; spin=true)

# Spinful basis ordering follows `kron(orbital, σ0)`:
#   index → (orbital, spin)   1↔(A,↑)  2↔(A,↓)  3↔(B,↑)  4↔(B,↓)
const A_up, A_dn, B_up, B_dn = 1, 2, 3, 4
const NN_offsets = [[0, 0], [-1, 0], [0, -1]]
const seed_channel = :chiral_plus  # :extended_s, :chiral_plus, or :chiral_minus

# Initial guess.
#
# Electron sector: ½ I — half-filled, spin-unpolarised, no off-diagonal density.
# A wholly-random ρ_init at amplitude O(1) confuses the BdG chempot bisection
# (the bands are unphysical until the SCF settles), and the iteration falls
# into a period-2 oscillation.
#
# Pairing sector: small NN singlet seed. For :chiral_plus the three C3-related
# bonds carry phases (1, ω, ω²), selecting the m = +1 chiral d+id irrep; use
# :chiral_minus for m = -1, or :extended_s for a real equal-bond seed.
#
# Pauli antisymmetry. The anomalous block satisfies
#   F[L][i,j] = −F[−L][j,i]   (from ⟨c_j(0) c_i(L)⟩ = −⟨c_i(L) c_j(0)⟩).
# For the on-site key (L = 0) this reduces to F[0][i,j] = −F[0][j,i] within
# the same block; failing to enforce it lets the seed drift into a spurious
# symmetric component that artificially breaks C3 across the three bonds on
# the first iteration.
ρ_el_init = zero(v)
ρ_el_init[TightBinding.zerokey(ρ_el_init)] += 0.5 * Matrix{ComplexF64}(I, 4, 4)

const Δ_seed_amp = 0.05
const ω_seed = cis(2π / 3)

Δ_init = zero(v)
for (a, δR) in enumerate(NN_offsets)
    phase = if seed_channel === :chiral_plus
        ω_seed^(a - 1)
    elseif seed_channel === :chiral_minus
        conj(ω_seed)^(a - 1)
    elseif seed_channel === :extended_s
        1 + 0im
    else
        error("Unknown seed_channel=$(seed_channel)")
    end
    amp = Δ_seed_amp * phase

    Δ_init[δR][B_up, A_dn] = +amp   # ⟨c_{A↓}(0) c_{B↑}(δR)⟩
    Δ_init[δR][B_dn, A_up] = -amp   # ⟨c_{A↑}(0) c_{B↓}(δR)⟩
    Δ_init[-δR][A_dn, B_up] = -amp  # Pauli partner
    Δ_init[-δR][A_up, B_dn] = +amp
end

ρ_init = BdGOperator(ρ_el_init, Δ_init)

# Full Hartree-Fock SCF. In fock-only mode the onsite U mostly penalises
# onsite anomalous density and both the real and chiral NN seeds flow to the
# normal fixed point. Retaining the density channel changes the self-consistent
# normal background and lets the chiral NN solution survive.
ρ_sol, ϵ_GS, mf, converged, residual = Meanfield.solvehartreefock(
    hbdg, v, ρ_init, filling;
    klin           = klin_scf,
    iterations     = 400,
    tol            = 1e-7,
    T              = T_scf,
    β              = 0.5,
    acceleration   = :anderson,
    anderson_depth = 5,
    verbose        = true,
)

@info "BdG SCF" converged residual ϵ_GS μ=mf.μ ϵ_H=mf.ϵH ϵ_F=mf.ϵF ϵ_pair=mf.ϵP
converged || error("BdG SCF did not converge; pairing analysis would be unreliable.")

# Quasiparticle bands. The tutorial pattern is:
#   hmf = Meanfield.hMF(HMF); Operators.addchemicalpotential!(hmf, -HMF.μ)
# before plotting. For BdG this shift must be applied to the full Nambu
# operator; `Superconductivity.electron` supplies the particle-sector weight
# used to color the particle-hole bands.
H_qp = Meanfield.hMF(mf)
Operators.addchemicalpotential!(H_qp, -mf.μ)
electron_weight = Superconductivity.electron(H_qp)
qp_bands = getbands(H_qp, ks_path, electron_weight)

p_qp = plot(qp_bands;
            ylabel         = "\$E_{\\mathrm{BdG}}/t\$",
            colorbar_title = "electron weight",
            markersize     = 2.0,
            markercolor    = :viridis,
            ylims          = (-0.75, 0.75),
            size           = (380, 300),
            title          = "BdG quasiparticles")
hline!(p_qp, [0.0]; lc=:black, lw=1, ls=:dash, label="")

# ---------------------------------------------------------------------------
# Step 3 — extract NN bond pairings and check C3 phase winding
# ---------------------------------------------------------------------------
println(":: Step 3 — read off NN pairings and check C3 winding")

# `getpairingsector` returns the upper-right Nambu block of the BdG density
# matrix as a `Hops`. Its (i, j) entry at offset δR is the Fourier component
# of the anomalous correlator ⟨c c⟩ (the pairing density) on the link
# B-orbital(δR) ↔ A-orbital(0).
ρΔ = Superconductivity.getpairingsector(ρ_sol)

# Each NN bond is identified by (δR_carrying_B, bond_vector_in_Cartesian).
NN_bonds = [
    (NN_offsets[1], (+1.0,  0.0       )),   # δ₁
    (NN_offsets[2], (-0.5, +sqrt(3)/2)),    # δ₂  = R(2π/3) δ₁
    (NN_offsets[3], (-0.5, -sqrt(3)/2)),    # δ₃  = R(4π/3) δ₁
]

# For each bond decompose the four spin amplitudes into singlet / triplet
# channels. With NN attractive V and no Hund-like coupling, the leading
# channel is singlet — but report all four to make any triplet admixture
# visible.
println()
println("NN bond pairing amplitudes from the converged ρ_BdG:")
println(rpad("  δR (cell)", 14), rpad("δ (Cart.)", 24),
        rpad("|Δ_singlet|", 14), rpad("arg/π", 11),
        rpad("|Δ_triplet0|", 16), "|Δ_↑↑|, |Δ_↓↓|")

singlets = ComplexF64[]
triplets = ComplexF64[]
for (δR, δvec) in NN_bonds
    blk = ρΔ[δR]
    Δ_uu = blk[B_up, A_up]                                  # equal-spin ↑↑
    Δ_dd = blk[B_dn, A_dn]                                  # equal-spin ↓↓
    Δ_s  = (blk[B_up, A_dn] - blk[B_dn, A_up]) / sqrt(2)    # singlet
    Δ_t0 = (blk[B_up, A_dn] + blk[B_dn, A_up]) / sqrt(2)    # triplet S_z=0
    push!(singlets, Δ_s)
    push!(triplets, Δ_t0)
    @printf("  (%2d,%2d)      (%+5.2f,%+5.2f)        %10.4e  %+7.3f      %10.4e      %.2e, %.2e\n",
            δR[1], δR[2], δvec[1], δvec[2],
            abs(Δ_s), angle(Δ_s)/π,
            abs(Δ_t0),
            abs(Δ_uu), abs(Δ_dd))
end

# Pick the dominant channel and analyse its C3 winding.
singlet_norm = sqrt(sum(abs2, singlets))
triplet_norm = sqrt(sum(abs2, triplets))
dominant, dom_label = singlet_norm >= triplet_norm ? (singlets, "singlet") :
                                                     (triplets, "triplet S_z=0")
dominant_norm = max(singlet_norm, triplet_norm)
println()
@printf("Dominant channel: %s  (‖Δ_s‖ = %.3e, ‖Δ_t0‖ = %.3e)\n",
        dom_label, singlet_norm, triplet_norm)

# C3 anisotropy: ratio of the largest to smallest |Δ| across the three bonds.
mags  = abs.(dominant)
ratio = maximum(mags) / max(minimum(mags), eps())
@printf("Bond magnitudes |Δ| = (%.3e, %.3e, %.3e)   max/min = %.2f\n",
        mags[1], mags[2], mags[3], ratio)

# Phase increments under successive C3 steps δ_a → δ_{a+1}; circular mean for
# robustness against per-bond phase jitter.
phase_steps = [angle(dominant[mod1(a + 1, 3)] / dominant[a]) for a in 1:3]
mean_step   = angle(sum(exp(1im * φ) for φ in phase_steps))   # ∈ (-π, π]
m_winding   = mean_step / (2π / 3)                            # m ≈ -1, 0, +1

@printf("Δarg per C3 step (deg): [%+7.2f, %+7.2f, %+7.2f]\n",
        rad2deg(phase_steps[1]), rad2deg(phase_steps[2]), rad2deg(phase_steps[3]))
@printf("⇒ C3 angular momentum m ≈ %+5.2f   (rounded: m = %+d)\n",
        m_winding, round(Int, m_winding))

if dominant_norm < 1e-5
    println("\nInterpretation: the anomalous density collapsed to the normal-state fixed")
    println("point. The phase winding is not meaningful when the pairing amplitude is")
    println("this small.")
elseif abs(round(Int, m_winding)) == 1
    println("\nInterpretation: m = ±1 chiral pairing — the singlet has a clean phase")
    println("winding of 2π/3 between successive C3-related bonds. With TR broken in the")
    println("normal state (Haldane), this is the d±id channel.")
elseif ratio > 2
    println("\nInterpretation: m ≈ 0 with strongly C3-anisotropic |Δ| — the SCF found a")
    println("nematic / Kekulé-like singlet. The lattice C3 symmetry is spontaneously")
    println("broken by the pairing pattern (one NN bond is weakened, the other two are")
    println("enhanced and equal — a residual C2 symmetry remains). This is NOT chiral")
    println("d±id; the chiral channel is higher-energy here or below the")
    println("finite-temperature/finite-grid pairing threshold for this setup.")
else
    println("\nInterpretation: m ≈ 0 with C3-symmetric |Δ| (max/min ≈ 1) — the SCF found")
    println("the trivial extended-s channel. The chiral d±id channel is degenerate or")
    println("not selected at this filling.")
end

# ---------------------------------------------------------------------------
# Plot the three Δ(δ_a) in the complex plane (right) and as a bar chart of |Δ|
# ---------------------------------------------------------------------------
maxmag = max(maximum(abs, dominant), eps())

p_phase = plot(; aspect_ratio = :equal,
               xlim = (-1.25 * maxmag, +1.25 * maxmag),
               ylim = (-1.25 * maxmag, +1.25 * maxmag),
               xlabel = "Re Δ", ylabel = "Im Δ",
               size  = (340, 340),
               title = @sprintf("Δ_%s on NN bonds  (m=%+d)",
                                replace(dom_label, " " => "_"),
                                round(Int, m_winding)),
               legend = :outerbottom, legend_columns = 3)

θ_ref = LinRange(0, 2π, 200)
plot!(p_phase, maxmag .* cos.(θ_ref), maxmag .* sin.(θ_ref);
      lc=:gray, ls=:dash, lw=1, label="")

bond_colors = [:tomato, :seagreen, :royalblue]
bond_labels = ["δ₁ (1,0)", "δ₂ (−1,0)", "δ₃ (0,−1)"]
for (a, Δ) in enumerate(dominant)
    plot!(p_phase, [0.0, real(Δ)], [0.0, imag(Δ)];
          lc=bond_colors[a], lw=3, label="")
    scatter!(p_phase, [real(Δ)], [imag(Δ)];
             mc=bond_colors[a], ms=8, markerstrokewidth=0,
             label=bond_labels[a])
end

# ---------------------------------------------------------------------------
# Save plots
# ---------------------------------------------------------------------------
mkpath("output")
savefig(p_bands, "output/haldane_chiral_sc_bands.pdf")
savefig(p_qp, "output/haldane_chiral_sc_quasiparticles.pdf")
savefig(p_phase, "output/haldane_chiral_sc_pairing.pdf")
savefig(plot(p_bands, p_qp, p_phase; layout=(1, 3), size=(1160, 340)),
        "output/haldane_chiral_sc.pdf")

println("\n:: Wrote output/haldane_chiral_sc{,_bands,_quasiparticles,_pairing}.pdf")

# =============================================================================
# Example (ii): charge–charge correlator inside a Wilson-line-quenched state.
#
# A spatial Wilson line W = χ†_{l1}[string]χ_{l2} creates a ±q charge pair joined by an
# electric-flux string — a "Wilson-line quench" of the vacuum, |ψ⟩ = W|vac⟩ (normalized).
# |ψ⟩ is NOT an energy eigenstate, so it evolves. We measure the CONNECTED, EQUAL-TIME
# current–current correlator inside that evolving state,
#
#     C(R, t) = ⟨ψ(t)| j¹(c+R) j¹(c−R) |ψ(t)⟩_c ,     |ψ(t)⟩ = e^{-iHt} W|vac⟩ ,
#
# for a set of detector separations R (bonds symmetric about the lattice centre c). This
# is an expectation value *within* the quench state — the correct object for a
# non-eigenstate (unlike the vacuum `correlator2pt`, which needs an eigenstate reference).
#
# Uses toolkit `wilson_line_quench` (→ act + normalize) and `detector_2pt_equal_time`
# (the exported operator product ChargeCurrent(bR)·ChargeCurrent(bL) + connected
# subtraction of ⟨j¹(bR)⟩⟨j¹(bL)⟩), evaluated on each TDVP snapshot.
#
# Run from the project root:
#   SMOKE=1 julia --project=. examples/charge_charge_wilson.jl
#   AG=0.2 MG=0.0 N=128 NL=16 T=12 RLIST=8,16,24 julia --project=. examples/charge_charge_wilson.jl
# =============================================================================
using Schwinger, SchwingerDetectors, JLD2, Printf

const SMOKE = get(ENV, "SMOKE", "0") == "1"
AG   = parse(Float64, get(ENV, "AG", "0.2"))
MG   = parse(Float64, get(ENV, "MG", "0.0"))
N    = parse(Int,     get(ENV, "N",  SMOKE ? "16" : "128"))
NL   = parse(Int,     get(ENV, "NL", SMOKE ? "4"  : "16"))     # Wilson-line length (sites)
T    = parse(Float64, get(ENV, "T",  SMOKE ? "1.0" : "12.0"))
NST  = parse(Int,     get(ENV, "NSTEPS", SMOKE ? "5" : "120"))
MAXB = parse(Int,     get(ENV, "MAXBOND", SMOKE ? "32" : "200"))
# detector separations R (bonds from centre): bond pair (c+R, c−R)
Rlist = [parse(Int, x) for x in split(get(ENV, "RLIST", SMOKE ? "2,4" : "8,16,24"), ",")]
OUT  = get(ENV, "OUT", joinpath(@__DIR__, "..", "..", "..", "data", "toolkit_examples"))
mkpath(OUT)

# 1. Vacuum, then the Wilson-line quench |ψ⟩ = W|vac⟩.
m   = build_model(; N = N, ag = AG, mg = MG)
c   = N ÷ 2
l1  = c - NL ÷ 2; l2 = l1 + NL
wl  = wilson_line_quench(m, l1, l2; gs_kwargs = (; energy_tol = 1e-9))   # normalized W|vac⟩
st0 = wl.state

# 2. Detector bond pairs, symmetric about the centre.
Rpairs = [(c + R, c - R) for R in Rlist]
@assert all(1 ≤ bL && bR ≤ N - 1 for (bR, bL) in Rpairs) "a detector bond is out of range"
@printf("wilson line on sites %d…%d;  ⟨j¹(c+R) j¹(c−R)⟩_c at R = %s\n", l1, l2, Rlist)

# 3. Connected equal-time ⟨j¹ j¹⟩ inside the evolving quench state, at every TDVP step.
cc(ψ) = (d = detector_2pt_equal_time(ψ, Rpairs; kind = :current, connected = true);
         [d[p] for p in Rpairs])
C0 = cc(st0)                                                            # t = 0 (the quench)
_, obs = evolve(st0, T; nsteps = NST, two_site = true, maxlinkdim = MAXB,
                observable = Dict("cc" => (ψ, t) -> cc(ψ)))
tvec = vcat(0.0, times(obs))
Cmat = vcat(permutedims(C0), density_map(obs, :cc))                    # (nt+1) × nR

jldsave(joinpath(OUT, "charge_charge_wilson$(SMOKE ? "_smoke" : "").jld2");
        N = N, ag = AG, mg = MG, l1 = l1, l2 = l2, Rlist = Rlist,
        t = tvec, C = Cmat)                                            # C[time, R-index]
@printf("C(R=%d): t=0 → %+.3e,   t=T → %+.3e\n", Rlist[1], Cmat[1, 1], Cmat[end, 1])
println("wrote $(joinpath(OUT, "charge_charge_wilson$(SMOKE ? "_smoke" : "").jld2"))")

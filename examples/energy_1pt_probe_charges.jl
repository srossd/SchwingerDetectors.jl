# =============================================================================
# Example (i): energy one-point function after inserting probe charges.
#
# Insert a static ±q probe-charge pair into the vacuum (a sudden "string" quench, done as
# a θ/2π step between the charges), then watch the ENERGY CURRENT 𝒥 = T⁰¹ radiate outward
# at two fixed detectors — i.e. the one-point function ⟨𝒥(x_D, t)⟩ of the energy detector.
#
# Uses the new built-ins: `string_quench` (→ Schwinger's `rehost`/`quench`), the
# `standard_densities` menu via `detector_1pt`, and `energy_densities(…; convention=:bond)`
# (the density that actually partners 𝒥 in ∂_t h_n = 𝒥_n − 𝒥_{n+1}).
#
# Run from the project root:
#   SMOKE=1 julia --project=examples examples/energy_1pt_probe_charges.jl
#   AG=0.2 MG=0.0 N=128 STRING_L=16 DET=24 T=14 julia --project=examples examples/energy_1pt_probe_charges.jl
# =============================================================================
using Schwinger, SchwingerDetectors, JLD2, Printf

const SMOKE = get(ENV, "SMOKE", "0") == "1"
AG   = parse(Float64, get(ENV, "AG", "0.2"))
MG   = parse(Float64, get(ENV, "MG", "0.0"))
N    = parse(Int,     get(ENV, "N",        SMOKE ? "16" : "128"))
Lstr = parse(Int,     get(ENV, "STRING_L", SMOKE ? "4"  : "16"))   # probe-charge separation (links)
DET  = parse(Int,     get(ENV, "DET",      SMOKE ? "3"  : "24"))   # detector offset from centre (sites)
T    = parse(Float64, get(ENV, "T",        SMOKE ? "1.0" : "12.0"))
NST  = parse(Int,     get(ENV, "NSTEPS",   SMOKE ? "5"  : "120"))
MAXB = parse(Int,     get(ENV, "MAXBOND",  SMOKE ? "32" : "200"))
OUT  = get(ENV, "OUT", joinpath(@__DIR__, "..", "..", "..", "data", "toolkit_examples"))
mkpath(OUT)

# 1. Prepare: vacuum, then a sudden ±q probe-charge string quench.
m = build_model(; N = N, ag = AG, mg = MG)                    # finite θ=0 lattice
q = string_quench(m; length_links = Lstr, gs_kwargs = (; energy_tol = 1e-9))
c = N ÷ 2
detsites = [c - DET, c + DET]                                 # one detector each side
@printf("probe charges on links %s;  energy detectors 𝒥 at sites %s\n", q.links, detsites)

# 2. Evolve under the string Hamiltonian and read the energy-current 1-pt at the detectors.
t, series = detector_1pt(q.state, T, detsites; kind = :energycurrent,
                         nsteps = NST, two_site = true, maxbond = MAXB)

# (also grab the bond-consistent energy-density spacetime map, the 𝒥 partner)
evolved, obs = evolve_and_measure(q.state, T; nsteps = NST, two_site = true, maxbond = MAXB,
                                  densities = [:energycurrent])
hbond = [energy_densities(evolved; convention = :bond)]       # single final snapshot, for illustration

jldsave(joinpath(OUT, "energy_1pt$(SMOKE ? "_smoke" : "").jld2");
        N = N, ag = AG, mg = MG, links = collect(q.links), detsites = detsites,
        t = t, J_left = series[detsites[1]], J_right = series[detsites[2]],
        Jmap = density_map(obs, :energycurrent))
@printf("⟨𝒥⟩ at right detector (site %d): t=0 → %.3e,   t=T → %.3e\n",
        detsites[2], series[detsites[2]][1], series[detsites[2]][end])
println("wrote $(joinpath(OUT, "energy_1pt$(SMOKE ? "_smoke" : "").jld2"))")

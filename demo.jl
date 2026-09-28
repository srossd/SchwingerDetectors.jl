# =============================================================================
# demo.jl — end-to-end smoke of the toolkit's three capabilities on a tiny lattice.
#
# Run from the project root (environment active):
#   SMOKE=1 julia --project=. Schwinger/toolkit/demo.jl        # tiny, ~seconds
#   AG=0.2 MG=0.0 N=64 julia --project=. Schwinger/toolkit/demo.jl
#
# It (1) builds a string quench, (2) evolves while measuring densities + reads a
# detector 1-pt time series, and (3) computes a detector 2-pt correlator — all with
# exported Schwinger.jl calls via SchwingerToolkit. Outputs go to data/toolkit_demo/.
# =============================================================================
using Schwinger, JLD2, Printf
include(joinpath(@__DIR__, "SchwingerToolkit.jl"))
using .SchwingerToolkit

const SMOKE = get(ENV, "SMOKE", "0") == "1"
const AG    = parse(Float64, get(ENV, "AG", "0.2"))
const MG    = parse(Float64, get(ENV, "MG", "0.0"))
const N     = parse(Int,     get(ENV, "N",  SMOKE ? "16" : "64"))
const T     = parse(Float64, get(ENV, "T",  SMOKE ? "1.0" : "6.0"))
const NST   = parse(Int,     get(ENV, "NSTEPS", SMOKE ? "5" : "60"))
const MAXB  = parse(Int,     get(ENV, "MAXBOND", SMOKE ? "32" : "128"))
const OUT   = get(ENV, "OUT", joinpath(@__DIR__, "..", "..", "data", "toolkit_demo"))
mkpath(OUT)

@printf("== toolkit demo: N=%d ag=%.3f m/g=%.3f T=%.1f nsteps=%d maxbond=%d ==\n",
        N, AG, MG, T, NST, MAXB); flush(stdout)

# ---- 1. state preparation: static ±q string quench ----
m  = build_model(; N = N, ag = AG, mg = MG)                      # finite θ=0 vacuum
Lstr = max(2, N ÷ 8)
q  = string_quench(m; length_links = Lstr, gs_kwargs = (; energy_tol = 1e-9))
c  = N ÷ 2
@printf("string quench: links %s, detectors at sites %d,%d\n", q.links, c - Lstr, c + Lstr)

# ---- 2. evolve while measuring densities, and read a detector 1-pt series ----
detsites = [c - Lstr, c + Lstr]
tvec, series = detector_1pt(q.state, T, detsites;
                            kind = :energycurrent, nsteps = NST,
                            two_site = true, maxbond = MAXB)
# also grab the full density spacetime maps for a couple of fields
evolved, obs = evolve_and_measure(q.state, T; nsteps = NST, two_site = true, maxbond = MAXB,
                                  densities = [:energy, :electricfield, :current, :energycurrent, :pseudoscalar])
jldsave(joinpath(OUT, "onept$(SMOKE ? "_smoke" : "").jld2");
        N = N, ag = AG, mg = MG, links = collect(q.links), detsites = detsites,
        times = tvec, det_energycurrent = Dict(string(s) => series[s] for s in detsites),
        energy = density_map(obs, :energy), electricfield = density_map(obs, :electricfield),
        current = density_map(obs, :current), energycurrent = density_map(obs, :energycurrent),
        pseudoscalar = density_map(obs, :pseudoscalar))
@printf("1-pt detector 𝒥 at site %d: first=%.3e last=%.3e\n",
        detsites[1], series[detsites[1]][1], series[detsites[1]][end])

# ---- 3. detector 2-pt correlator ⟨j¹(bR,t) j¹(bL,0)⟩ on the vacuum ----
lat  = m.lat
bL   = c - Lstr
bR   = c + Lstr
tgrid = collect(0.0:(T / NST):T)
tc, C = detector_2pt(q.gs,
                     ChargeCurrent(lat, bL; backend = :MPSKit),   # source j¹ at bL
                     ChargeCurrent(lat, bR; backend = :MPSKit),   # detector j¹ at bR
                     tgrid; connected = true, two_site = true, maxlinkdim = MAXB)
jldsave(joinpath(OUT, "twopt$(SMOKE ? "_smoke" : "").jld2");
        N = N, ag = AG, mg = MG, bL = bL, bR = bR, t = tc, re = real.(C), im = imag.(C))
@printf("2-pt ⟨j¹(%d,t) j¹(%d,0)⟩_c: C(0)=%.3e  |C(T)|=%.3e\n", bR, bL, real(C[1]), abs(C[end]))

println("wrote outputs to $OUT")
println("=== demo done ===")

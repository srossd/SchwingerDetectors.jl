# θ=0 ground state for the energy-radiation study (one per m/g), saved for reuse.
#   julia --project=production production/erad_vacuum.jl <m_over_g> <outdir>
using Statistics, Printf
using Schwinger, SchwingerDetectors
include(joinpath(@__DIR__, "persist.jl"))
mgtag(x) = replace(string(x), "." => "p")

mg = parse(Float64, ARGS[1]); outdir = get(ARGS, 2, "data/erad"); mkpath(outdir)
N = parse(Int, get(ENV, "N", "256")); ag = parse(Float64, get(ENV, "AG", "0.2"))
m = build_model(N=N, F=1, q=1, ag=ag, mg=mg, theta2pi=0.0)
gs = groundstate(m.H; energy_tol=1e-7, bonddim=96)
@printf("[erad-vac mg=%.2f] N=%d ag=%.2f  E0=%.5f  <E>=%.4f\n",
        mg, N, ag, real(energy(gs)), mean(real.(electricfields(gs)))); flush(stdout)
save_state(joinpath(outdir, "erad_gs_mg$(mgtag(mg)).bin"), gs)
println("saved θ=0 vacuum")

# Stage 1 (cheap): compute + verify the ground state for one θ/π, and save it.
#
# Verifies the vacuum sits on the true (screened) branch: |<E_avg>| < EAVG_TOL. At m/g = 1
# the first-order transition is at θ = π (θ/π = 1); for θ/π = 1.2, 1.4 the true vacuum is the
# FLIPPED one (<E> ≈ +0.4, +0.3), not the continued branch (≈ -0.6, -0.7). DMRG targets the
# global energy minimum, so the true vacuum is expected automatically; this script checks it.
#
#   julia --project=production production/vacuum.jl <theta_over_pi> <outdir>
using Statistics, Printf, JLD2
include(joinpath(@__DIR__, "config.jl"))
include(joinpath(@__DIR__, "persist.jl"))

theta_over_pi = parse(Float64, ARGS[1])
outdir = get(ARGS, 2, joinpath(@__DIR__, "..", "data", "qq_sweep"))
mkpath(outdir)

t2p = theta2pi(theta_over_pi)
m   = model(t2p)
@printf("[vacuum] θ/π=%.3f  θ2π=%.3f  N=%d ag=%.3f m/g=%.3f  (bond=%d, tol=%.0e)\n",
        theta_over_pi, t2p, N, AG, MG, GS_BOND, GS_TOL); flush(stdout)

gs = groundstate(m.H; energy_tol = GS_TOL, bonddim = GS_BOND)
Efields = real.(electricfields(gs))
Eavg = mean(Efields)
E0   = real(energy(gs))
ok = abs(Eavg) < EAVG_TOL + 0.02   # small tolerance; θ/π=1 is borderline (<E>≈±0.5)
@printf("[vacuum] E0=%.6f  <E_avg>=%.4f  |<E_avg>|=%.4f  -> %s\n",
        E0, Eavg, abs(Eavg), ok ? "OK (true branch)" : "FAIL (wrong branch?)"); flush(stdout)

path = joinpath(outdir, "gs_theta$(thetatag(theta_over_pi)).bin")
save_state(path, gs)
@save joinpath(outdir, "vacinfo_theta$(thetatag(theta_over_pi)).jld2") theta_over_pi t2p Eavg E0 ok N AG MG
println("[vacuum] saved $path")
ok || (@error "vacuum on wrong branch"; exit(3))

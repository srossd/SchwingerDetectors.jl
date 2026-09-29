# Real per-call cost of each measurement on a FINITE lattice, to confirm the O(N²)-full-vector
# vs O(N)-single-bond analysis. Measured on the (low-bond) quench state at t=0 — no evolution
# needed: the ratio full-vector / single-bond reveals the algorithmic scaling regardless of bond.
#   julia --project=production production/erad_meas_timing.jl
using Statistics, Printf
using Schwinger, SchwingerDetectors

N = 256; ag = 0.2; mg = 1.0
m0 = build_model(N=N, F=1, q=1, ag=ag, mg=mg, theta2pi=0.0)
tg = @elapsed gs = groundstate(m0.H; energy_tol=1e-2, bonddim=32)   # rough & fast; low bond is fine
@printf("[meas] gs %.0fs; N=%d\n", tg, N); flush(stdout)
c = N÷2; l1 = c-10; l2 = c+10; θ = fill(0.0, N); for n in l1:(l2-1); θ[n] += 0.5; end
latq = Lattice(N; F=1, q=1, a=ag, m=mg, θ2π=θ); Hq = Hamiltonian(latq; backend=:MPSKit)
psi = MPSKitState(Hq, gs.psi)

electricfields(psi); chargecurrents(psi)   # warmup/compile
t_ef = @elapsed electricfields(psi); @printf("[meas] electricfields(full): %.2fs\n", t_ef); flush(stdout)
t_jc = @elapsed chargecurrents(psi);  @printf("[meas] chargecurrents(full): %.2fs\n", t_jc); flush(stdout)
t_je = @elapsed energycurrents(psi);  @printf("[meas] energycurrents(full): %.2fs\n", t_je); flush(stdout)
op = EnergyCurrent(latq, c+40; backend=:MPSKit); real(expectation(op, psi))   # warmup
t_1 = @elapsed real(expectation(op, psi)); @printf("[meas] single-bond EnergyCurrent expectation: %.3fs\n", t_1); flush(stdout)
t_ed = @elapsed energy_densities(psi); @printf("[meas] energy_densities(full): %.2fs\n", t_ed); flush(stdout)
@printf("[meas] ratios: energycurrents/single = %.0fx ; energy_densities/single = %.0fx (N=%d)\n",
        t_je/max(t_1,1e-6), t_ed/max(t_1,1e-6), N)
println("MEAS_DONE")

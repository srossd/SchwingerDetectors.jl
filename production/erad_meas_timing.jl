# Real per-call cost of each measurement on a FINITE lattice (moderate bond), to confirm the
# O(N²)-full-vector vs O(N)-single-bond analysis. N=256, θ=π quench, a few steps for bond growth.
#   julia --project=production production/erad_meas_timing.jl
using Statistics, Printf
using Schwinger, SchwingerDetectors

N = 256; ag = 0.2; mg = 1.0
m0 = build_model(N=N, F=1, q=1, ag=ag, mg=mg, theta2pi=0.0)
gs = groundstate(m0.H; energy_tol=1e-3, bonddim=64)
c = N÷2; l1 = c-10; l2 = c+10; θ = fill(0.0, N); for n in l1:(l2-1); θ[n] += 0.5; end
latq = Lattice(N; F=1, q=1, a=ag, m=mg, θ2π=θ); Hq = Hamiltonian(latq; backend=:MPSKit)
psi = MPSKitState(Hq, gs.psi)
for k in 1:4; global psi = evolve(psi, 0.2; nsteps=1, two_site=true, maxlinkdim=256)[1]; end   # grow some bond
@printf("[meas] N=%d, after 4 θ=π steps\n", N); flush(stdout)

t_ef = @elapsed ef = electricfields(psi)
t_jc = @elapsed jc = chargecurrents(psi)
t_je = @elapsed je = energycurrents(psi)
t_ed = @elapsed ed = energy_densities(psi)
op = EnergyCurrent(latq, c+40; backend=:MPSKit)
t_1 = @elapsed v = real(expectation(op, psi))
@printf("[meas] electricfields(full) %.2fs | chargecurrents(full) %.2fs | energycurrents(full) %.2fs | energy_densities(full) %.2fs\n",
        t_ef, t_jc, t_je, t_ed)
@printf("[meas] SINGLE-bond EnergyCurrent expectation: %.3fs  (ratio full/single = %.0fx)\n", t_1, t_je/max(t_1,1e-6))
println("MEAS_DONE")

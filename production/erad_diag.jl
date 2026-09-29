# Isolate why the erad runs are pathologically slow: scalar-θ vs θ-vector (string) TDVP step
# cost, Hamiltonian build time, and the cost of the energy_densities/energycurrents/electricfields
# measurements — all at N=256, m/g=1 (fast gs).
#   julia --project=production production/erad_diag.jl
using Statistics, Printf
using Schwinger, SchwingerDetectors

N = 256; ag = 0.2; mg = 1.0; maxb = 256
step(s) = evolve(s, 0.2; nsteps=1, two_site=true, maxlinkdim=maxb)[1]

m0 = build_model(N=N, F=1, q=1, ag=ag, mg=mg, theta2pi=0.0)
tgs = @elapsed gs = groundstate(m0.H; energy_tol=1e-4, bonddim=48)
@printf("[diag] rough gs (m/g=1, N=%d): %.0fs\n", N, tgs); flush(stdout)

# --- scalar-θ TDVP steps (like bench512/QQ) ---
ss = step(gs)                              # warmup/compile
t_scalar = @elapsed for _ in 1:3; global ss = step(ss); end
@printf("[diag] scalar-θ  step: %.1f s/step\n", t_scalar/3); flush(stdout)

# --- θ-vector (string) Hamiltonian build + TDVP steps ---
c = N ÷ 2; l1 = c-10; l2 = c+10
θ = fill(0.0, N); for n in l1:(l2-1); θ[n] += 0.5; end
tH = @elapsed Hq = Hamiltonian(Lattice(N; F=1, q=1, a=ag, m=mg, θ2π=θ); backend=:MPSKit)
@printf("[diag] θ-vector Hamiltonian build: %.1f s\n", tH); flush(stdout)
sv = MPSKitState(Hq, gs.psi)
sv = step(sv)                              # warmup/compile
t_vec = @elapsed for _ in 1:3; global sv = step(sv); end
@printf("[diag] θ-vector step: %.1f s/step\n", t_vec/3); flush(stdout)

# --- measurement costs (called every STRIDE in erad_run) ---
t_ed = @elapsed energy_densities(sv); t_ec = @elapsed energycurrents(sv); t_ef = @elapsed electricfields(sv)
@printf("[diag] measurements: energy_densities %.1fs, energycurrents %.1fs, electricfields %.1fs\n", t_ed, t_ec, t_ef)
println("DIAG_DONE")

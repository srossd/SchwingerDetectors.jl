# Timing/feasibility probe on the real 1024-site problem, so dt/maxbond/walltime/task
# granularity can be chosen with evidence (not guessed). Times the groundstate solve, a
# savestate/loadstate round-trip, and a handful of streaming-correlator steps for 1 and 3
# detector pairs, then extrapolates to the full run.
#
#   julia --project=production production/benchmark.jl [theta_over_pi]
using Statistics, Printf
include(joinpath(@__DIR__, "config.jl"))

theta_over_pi = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 1.2
t2p = theta2pi(theta_over_pi)
m = model(t2p)
@printf("[bench] N=%d ag=%.3f m/g=%.3f θ/π=%.2f  MAXBOND=%d SUBSTEPS=%d DT=%.3f\n",
        N, AG, MG, theta_over_pi, MAXBOND, SUBSTEP, DT); flush(stdout)

tgs = @elapsed gs = groundstate(m.H; energy_tol = GS_TOL, bonddim = GS_BOND)
Eavg = mean(real.(electricfields(gs)))
@printf("[bench] groundstate: %.1f s   <E_avg>=%.4f\n", tgs, Eavg); flush(stdout)

# savestate/loadstate round-trip
tmp = tempname() * ".jld2"
savestate(tmp, gs)
gs2 = loadstate(tmp, m.H)
@printf("[bench] savestate/loadstate ok: |<E>-<E'>|=%.2e\n",
        abs(Eavg - mean(real.(electricfields(gs2))))); flush(stdout)

l1, l2 = wilson_endpoints()
wl = wilson_line_quench(m, l1, l2; gs = gs)
psi = wl.state
@printf("[bench] wilson line sites %d…%d (gx=%.1f, %d sites)\n", l1, l2, GX_WL, l2 - l1); flush(stdout)

K = 4
onepair = [(CENTER - RLIST[1], CENTER + RLIST[1])]
threepair = [(CENTER - R, CENTER + R) for R in RLIST]

t1 = @elapsed charge_transfer_correlator(psi, K*DT, onepair; nsteps=K, substeps=SUBSTEP, two_site=true, maxbond=MAXBOND)
t3 = @elapsed charge_transfer_correlator(psi, K*DT, threepair; nsteps=K, substeps=SUBSTEP, two_site=true, maxbond=MAXBOND)
Mfull = round(Int, T / DT)
@printf("[bench] %d steps: 1-pair %.1f s (%.2f s/step), 3-pair %.1f s (%.2f s/step)\n",
        K, t1, t1/K, t3, t3/K); flush(stdout)
@printf("[bench] EXTRAPOLATION to T=%.0f (M=%d): 3-pair single-task ≈ %.1f min = %.2f h\n",
        T, Mfull, t3/K*Mfull/60, t3/K*Mfull/3600); flush(stdout)
@printf("[bench] per (θ,R) task (1 pair) ≈ %.2f h\n", t1/K*Mfull/3600); flush(stdout)
println("[bench] done")

# Real cost of the ag=0.2, N=512 (same physical box) plan: per-step TDVP2 cost as entanglement
# grows, memory, and rough-gs time. Uses the folded θ2π=-0.4 (θ/π=1.2) true vacuum.
#   julia --project=production production/bench512.jl
using Statistics, Printf
using Schwinger, SchwingerDetectors, MPSKit

peakGB() = parse(Int, match(r"VmHWM:\s+(\d+)", read("/proc/self/status", String)).captures[1]) / 1024^2
maxent(s) = try; maximum(real.(entanglements(s))); catch; -1.0; end

N = 512; ag = 0.2
m = build_model(; N=N, ag=ag, mg=1.0, theta2pi=-0.4)
tgs = @elapsed g = groundstate(m.H; energy_tol=1e-3, bonddim=48)   # rough, for timing only
@printf("[b512] rough gs: %.0fs  <E>=%.4f  maxS=%.3f  peakRSS=%.2f GB\n",
        tgs, mean(real.(electricfields(g))), maxent(g), peakGB()); flush(stdout)

c = N ÷ 2
psi = wilson_line_quench(m, c-5, c+5; gs=g).state    # gx=2 -> 10 sites
dt = 0.2; maxb = 256
psi = evolve(psi, dt; nsteps=1, two_site=true, maxlinkdim=maxb)[1]   # warmup/compile
tsum = 0.0
for k in 1:16
    t = @elapsed global psi = evolve(psi, dt; nsteps=1, two_site=true, maxlinkdim=maxb)[1]
    global tsum += t
    @printf("[b512] step %2d (t=%.1f): %.1f s/step  maxS=%.3f  peakRSS=%.2f GB\n",
            k, k*dt, t, maxent(psi), peakGB()); flush(stdout)
end
M = round(Int, 90/dt)
@printf("[b512] avg %.1f s/step over 16 steps; naive full-run (5 states, M=%d) = %.1f h/task\n",
        tsum/16, M, 5*(tsum/16)*M/3600)
println("B512_DONE")

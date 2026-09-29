# Focused diagnostic: WHY is gs slow, and what does ONE TDVP step actually cost (scaling in N)?
# Decides feasibility of the T=90 evolution at large N.
#   julia --project=production production/diag2.jl
using Statistics, Printf
using Schwinger, SchwingerDetectors, MPSKit

peakGB() = parse(Int, match(r"VmHWM:\s+(\d+)", read("/proc/self/status", String)).captures[1]) / 1024^2

function bonds_of(state)
    ψ = state.psi
    try; return maximum(MPSKit.bond_dimensions(ψ)); catch; end
    try; return maximum(i -> dim(right_virtualspace(ψ, i)), 1:length(ψ)); catch; end
    return -1
end

# --- gs verbose at N=128: sweeps, time, <E>, MPS length (gauge sites?) ---
m128 = build_model(; N=128, ag=0.1, mg=1.0, theta2pi=0.6)
t = @elapsed gs = groundstate(m128.H; energy_tol=1e-8, bonddim=64, verbose=true)
@printf("[gs] N=128 mpslen=%d: %.1fs  <E>=%.4f  stateMaxBond=%d  peakRSS=%.2f GB\n",
        length(gs.psi), t, mean(real.(electricfields(gs))), bonds_of(gs), peakGB()); flush(stdout)

# --- single TDVP2 evolve-step cost vs N ---
for Ns in (128, 256, 512)
    m  = build_model(; N=Ns, ag=0.1, mg=1.0, theta2pi=0.6)
    g  = groundstate(m.H; energy_tol=1e-6, bonddim=32)
    c  = Ns ÷ 2
    ψ  = wilson_line_quench(m, c-10, c+10; gs=g).state
    evolve(ψ, 0.1; nsteps=1, two_site=true, maxlinkdim=128)          # warmup/compile
    ts = @elapsed for _ in 1:3; global ψ = evolve(ψ, 0.1; nsteps=1, two_site=true, maxlinkdim=128)[1]; end
    @printf("[step] N=%d: %.2f s/TDVP2-step (maxbond=128)  peakRSS=%.2f GB\n", Ns, ts/3, peakGB()); flush(stdout)
end
println("DIAG2_DONE")

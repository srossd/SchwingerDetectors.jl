# Isolate the TDVP-step cost (the number that decides evolution feasibility). Uses a rough,
# cheap state so gs cost doesn't dominate; timing depends on maxlinkdim, not gs quality.
#   julia --project=production production/diag3.jl
using Statistics, Printf
using Schwinger, SchwingerDetectors, MPSKit

peakGB() = parse(Int, match(r"VmHWM:\s+(\d+)", read("/proc/self/status", String)).captures[1]) / 1024^2

function time_steps(Ns, maxb)
    m = build_model(; N=Ns, ag=0.1, mg=1.0, theta2pi=0.6)
    tgs = @elapsed g = groundstate(m.H; energy_tol=1e-2, bonddim=16)   # deliberately rough/cheap
    c = Ns ÷ 2
    ψ = wilson_line_quench(m, c-10, c+10; gs=g).state
    ψ = evolve(ψ, 0.1; nsteps=1, two_site=true, maxlinkdim=maxb)[1]     # warmup/compile
    ts = @elapsed for _ in 1:3
        ψ = evolve(ψ, 0.1; nsteps=1, two_site=true, maxlinkdim=maxb)[1]
    end
    @printf("[step] N=%4d maxbond=%d: rough-gs %.0fs;  %.2f s/TDVP2-step;  peakRSS=%.2f GB\n",
            Ns, maxb, tgs, ts/3, peakGB()); flush(stdout)
end

for Ns in (128, 256, 512)
    time_steps(Ns, 128)
end
println("DIAG3_DONE")

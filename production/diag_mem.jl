# Memory diagnostic: find what actually drives RAM in the groundstate solve.
# Measures peak resident memory (VmHWM) and the bond dim the state actually needs, vs N and
# the requested groundstate bond cap. m/g=1 is gapped, so the true state should need small bond.
#
#   julia --project=production production/diag_mem.jl
using Statistics, Printf
using Schwinger, SchwingerDetectors

peakGB() = parse(Int, match(r"VmHWM:\s+(\d+)", read("/proc/self/status", String)).captures[1]) / 1024^2  # KB->GB

function maxbond_of(state)
    ψ = state.psi
    try
        return maximum(i -> dim(MPSKit.left_virtualspace(ψ, i)), 1:length(ψ))
    catch
        return -1
    end
end

for (Nv, B) in ((128, 50), (128, 200), (256, 100), (512, 100), (512, 50))
    lat = Lattice(Nv; F=1, q=1, a=0.1, m=1.0, θ2π=0.6)
    H = Hamiltonian(lat; backend=:MPSKit)
    t = @elapsed gs = groundstate(H; energy_tol=1e-6, bonddim=B)
    Eavg = mean(real.(electricfields(gs)))
    @printf("N=%4d bondcap=%3d : %.1fs  peakRSS=%.2f GB  stateMaxBond=%d  <E>=%.4f\n",
            Nv, B, t, peakGB(), maxbond_of(gs), Eavg); flush(stdout)
    GC.gc()
end
println("DIAG_DONE peak=", round(peakGB(); digits=2), " GB")

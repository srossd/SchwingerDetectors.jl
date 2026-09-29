# Test the vacuum-branch fix: θ2π = +0.6 and θ2π = -0.4 are the SAME physical θ (θ periodic 2π),
# but DMRG's random L≈0 init should land on the false branch for +0.6 and the true (screened)
# branch for -0.4. If so, folding θ2π into [-0.5,0.5] selects the true vacuum (|<E>|<0.5).
#   julia --project=production production/diag_vac.jl
using Statistics, Printf
using Schwinger, SchwingerDetectors
Nv = 96
for t2p in (0.6, -0.4, 0.7, -0.3)
    m = build_model(; N=Nv, ag=0.1, mg=1.0, theta2pi=t2p)
    t = @elapsed gs = groundstate(m.H; energy_tol=1e-6, bonddim=48)
    E = real.(electricfields(gs))
    @printf("θ2π=%+.2f (N=%d): %.0fs  <E>=%+.4f  |<E>|=%.4f  E0=%.5f  -> %s\n",
            t2p, Nv, t, mean(E), abs(mean(E)), real(energy(gs)),
            abs(mean(E))<0.5 ? "TRUE branch" : "false branch"); flush(stdout)
end
println("DIAGVAC_DONE")

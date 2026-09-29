# maxbond convergence check for the ±1/2 (θ=π) energy-radiation quench, worst case (m/g=0,
# highest entanglement). Evolve the same quench at several maxbonds and compare the radiated
# energy E_rad = ∫𝒥 dt at a fixed detector. Smaller N=192 for speed; the bond needed is set by
# the (local) entanglement near the string, ~N-independent.
#   julia --project=production production/erad_maxbond_check.jl
using Statistics, Printf
using Schwinger, SchwingerDetectors

N = 192; ag = 0.2; mg = 0.0; L = 15; T = 12.0; dt = 0.2; dtheta = 0.5
M = round(Int, T/dt); c = N ÷ 2
trapz(t, y) = sum((@view(y[2:end]) .+ @view(y[1:end-1]))/2 .* diff(t))

m0 = build_model(N=N, F=1, q=1, ag=ag, mg=mg, theta2pi=0.0)
tg = @elapsed gs = groundstate(m0.H; energy_tol=1e-4, bonddim=64)
@printf("[bcheck] rough massless gs (N=%d): %.0fs  <E>=%.4f\n", N, tg, mean(real.(electricfields(gs)))); flush(stdout)

l1 = c - L÷2; l2 = l1 + L; θ = fill(0.0, N); for n in l1:(l2-1); θ[n] += dtheta; end
latq = Lattice(N; F=1, q=1, a=ag, m=mg, θ2π=θ); Hq = Hamiltonian(latq; backend=:MPSKit)
Ddet = c + 40                                   # detector 40 sites (phys 8) from centre, outside string (±7)
jEop = EnergyCurrent(latq, Ddet; backend=:MPSKit)

for maxb in (96, 160, 256)
    psi = MPSKitState(Hq, gs.psi)
    ts = [0.0]; js = [real(expectation(jEop, psi))]
    tw = @elapsed for k in 1:M
        psi = evolve(psi, dt; nsteps=1, two_site=true, maxlinkdim=maxb)[1]
        push!(ts, k*dt); push!(js, real(expectation(jEop, psi)))
    end
    @printf("[bcheck] maxbond=%3d: E_rad(T=%.0f, x_D=phys8)=%.5e   (%.1f s/step)\n",
            maxb, T, trapz(ts, js), tw/M); flush(stdout)
end
println("BCHECK_DONE")

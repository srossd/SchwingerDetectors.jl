# Energy radiated after inserting a ±1/2 static probe-charge pair at separation L.
#
# Prepare the θ=0 vacuum, suddenly impose a θ2π = +dtheta (=1/2 for ±1/2 charges) step over L
# links about the centre (the Gauss-law field of a ±1/2 defect pair), and evolve under the
# stepped Hamiltonian. Record the 1-pt energy-current correlator ⟨𝒥(x,t)⟩ (energy flux) and the
# energy-density / electric-field spacetime maps. The radiated energy past a detector x_D is
# E_rad = ∫₀ᵀ ⟨𝒥(x_D,t)⟩ dt. Sweeping L tests confinement: E should saturate in L at m/g=0
# (screening, no linear term) but grow ~linearly at m/g=1 (½-charge string tension at θ=π).
#
# NB (cluster Schwinger = GitHub main): rehost() and standard_densities() are absent, so we
# "rehost" by hand as MPSKitState(H_step, gs.psi) and evolve with a manual observable using the
# vectorized energycurrents / energy_densities / electricfields.
#
#   julia --project=production production/erad_run.jl <m_over_g> <L_sites> <outdir>
using Statistics, Printf, JLD2
using Schwinger, SchwingerDetectors
include(joinpath(@__DIR__, "persist.jl"))
mgtag(x) = replace(string(x), "." => "p")

mg = parse(Float64, ARGS[1]); L = parse(Int, ARGS[2]); outdir = get(ARGS, 3, "data/erad"); mkpath(outdir)
N = parse(Int, get(ENV, "N", "256")); ag = parse(Float64, get(ENV, "AG", "0.2"))
T = parse(Float64, get(ENV, "T", "40.0")); dt = parse(Float64, get(ENV, "DT", "0.2"))
maxbond = parse(Int, get(ENV, "MAXBOND", "256")); dtheta = parse(Float64, get(ENV, "DTHETA", "0.5"))
M = round(Int, T / dt); c = N ÷ 2
outfile = joinpath(outdir, "erad_mg$(mgtag(mg))_L$(L).jld2")
isfile(outfile) && (println("already done: $outfile"); exit(0))

# θ=0 vacuum (load if present, else compute)
m0 = build_model(N=N, F=1, q=1, ag=ag, mg=mg, theta2pi=0.0)
gsf = joinpath(outdir, "erad_gs_mg$(mgtag(mg)).bin")
gs = isfile(gsf) ? load_state(gsf, m0.H) : (g = groundstate(m0.H; energy_tol=1e-7, bonddim=96); save_state(gsf, g); g)
@printf("[erad mg=%.2f L=%d] θ=0 vacuum <E>=%.4f\n", mg, L, mean(real.(electricfields(gs)))); flush(stdout)

# sudden ±dtheta probe-charge string over L links about the centre; rehost by hand
l1 = c - L ÷ 2; l2 = l1 + L
θ = fill(0.0, N); for n in l1:(l2 - 1); θ[n] += dtheta; end
latq = Lattice(N; F=1, q=1, a=ag, m=mg, θ2π=θ)
Hq = Hamiltonian(latq; backend=:MPSKit)
qstate = MPSKitState(Hq, gs.psi)               # = rehost(gs, Hq)
@printf("[erad mg=%.2f L=%d] string links %d…%d (gx=%.2f), dtheta=%.2f; T=%d M=%d dt=%.2f\n",
        mg, L, l1, l2, L*ag, dtheta, Int(T), M, dt); flush(stdout)

erow(s) = vec(real.(energycurrents(s))); hrow(s) = vec(real.(energy_densities(s))); efrow(s) = vec(real.(electricfields(s)))
hvac = hrow(gs)                                  # θ=0 vacuum energy-density reference
obs_dict = Dict("jE" => (ψ,t) -> erow(ψ), "h" => (ψ,t) -> hrow(ψ), "ef" => (ψ,t) -> efrow(ψ))
trun = @elapsed (evolved, obs) = evolve(qstate, T; nsteps=M, two_site=true, maxlinkdim=maxbond, observable=obs_dict)
@printf("[erad mg=%.2f L=%d] evolved in %.1f min\n", mg, L, trun/60); flush(stdout)

tg = times(obs); Jmap = density_map(obs, "jE"); Hmap = density_map(obs, "h"); EFmap = density_map(obs, "ef")
jldsave(outfile; mg=mg, L=L, ag=ag, N=N, c=c, T=T, dt=dt, dtheta=dtheta, links=[l1,l2],
        t=tg, Jmap=Jmap, Hmap=Hmap, EFmap=EFmap, hvac=hvac)
@printf("[erad mg=%.2f L=%d] injected E (Σ h−h_vac, t=0) = %.4f  -> %s\n",
        mg, L, sum(Hmap[1,:] .- hvac), outfile)

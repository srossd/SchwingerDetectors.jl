# Energy radiated after inserting a ±1/2 static probe-charge pair at separation L.
#
# Prepare the θ=0 vacuum, suddenly impose a θ2π=+dtheta (=1/2 for ±1/2 charges) step over L
# links about the centre, evolve under the stepped Hamiltonian, and record the 1-pt
# energy-current correlator ⟨𝒥(x,t)⟩ (energy flux) + energy-density / electric-field maps.
# Radiated energy past a detector x_D is E_rad = ∫₀ᵀ ⟨𝒥(x_D,t)⟩ dt.
#
# Explicit step loop with FLUSHED periodic progress + checkpoint/resume + wall budget, so
# progress is visible and a long (m/g=0, massless) run can't be silently lost.
# NB (cluster Schwinger = main): rehost()/standard_densities() absent -> rehost by hand as
# MPSKitState(H_step, gs.psi) and use vectorized energycurrents/energy_densities/electricfields.
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
STRIDE = parse(Int, get(ENV, "STRIDE", "2"))          # record maps every STRIDE steps
PRINT_EVERY = parse(Int, get(ENV, "PRINT_EVERY", "5"))
CKPT_EVERY = parse(Int, get(ENV, "CKPT_EVERY", "25"))
WALL = parse(Float64, get(ENV, "WALL_SECONDS", "61200"))   # ~17h
M = round(Int, T / dt); c = N ÷ 2
outfile = joinpath(outdir, "erad_mg$(mgtag(mg))_L$(L).jld2")
isfile(outfile) && (println("already done: $outfile"); exit(0))
ckd = joinpath(outdir, "erad_ckpt_mg$(mgtag(mg))_L$(L)"); mkpath(ckd)

# θ=0 vacuum (load if present, else compute)
m0 = build_model(N=N, F=1, q=1, ag=ag, mg=mg, theta2pi=0.0)
gsf = joinpath(outdir, "erad_gs_mg$(mgtag(mg)).bin")
gs = isfile(gsf) ? load_state(gsf, m0.H) : (g = groundstate(m0.H; energy_tol=1e-7, bonddim=96); save_state(gsf, g); g)

# ±dtheta probe-charge string over L links; rehost by hand under the stepped Hamiltonian
l1 = c - L ÷ 2; l2 = l1 + L
θ = fill(0.0, N); for n in l1:(l2 - 1); θ[n] += dtheta; end
Hq = Hamiltonian(Lattice(N; F=1, q=1, a=ag, m=mg, θ2π=θ); backend=:MPSKit)
erow(s) = vec(real.(energycurrents(s))); hrow(s) = vec(real.(energy_densities(s))); efrow(s) = vec(real.(electricfields(s)))
hvac = hrow(gs)

metaf = joinpath(ckd, "meta.jld2")
if isfile(metaf)
    d = load(metaf); k = d["k"]::Int
    heat_t = d["heat_t"]; Jc = d["Jc"]; Hd = d["Hd"]; Ef = d["Ef"]
    psi = load_state(joinpath(ckd, "psi.bin"), Hq)
    @printf("[erad mg=%.2f L=%d] RESUME at step %d/%d\n", mg, L, k, M); flush(stdout)
else
    psi = MPSKitState(Hq, gs.psi)                 # = rehost(gs, Hq)
    heat_t = [0.0]; Jc = [erow(psi)]; Hd = [hrow(psi)]; Ef = [efrow(psi)]; k = 0
    @printf("[erad mg=%.2f L=%d] FRESH: θ=0 <E>=%.4f; string %d…%d (gx=%.2f) dθ=%.2f; N=%d T=%d M=%d dt=%.2f maxbond=%d\n",
            mg, L, mean(real.(electricfields(gs))), l1, l2, L*ag, dtheta, N, Int(T), M, dt, maxbond); flush(stdout)
end

function checkpoint(k)
    save_state(joinpath(ckd, "psi.bin"), psi)
    tmp = metaf * ".tmp"; jldsave(tmp; k=k, heat_t=heat_t, Jc=Jc, Hd=Hd, Ef=Ef); mv(tmp, metaf; force=true)
end

t0 = time()
while k < M
    global k += 1
    global psi = evolve(psi, dt; nsteps=1, two_site=true, maxlinkdim=maxbond)[1]
    if k % STRIDE == 0
        push!(heat_t, k*dt); push!(Jc, erow(psi)); push!(Hd, hrow(psi)); push!(Ef, efrow(psi))
    end
    if k % PRINT_EVERY == 0
        el = time() - t0
        @printf("[erad mg=%.2f L=%d] step %d/%d  elapsed=%.2fh  proj-total=%.1fh\n",
                mg, L, k, M, el/3600, el/k*M/3600); flush(stdout)
    end
    if k % CKPT_EVERY == 0 || (time()-t0) > WALL
        checkpoint(k)
        if (time()-t0) > WALL && k < M
            @printf("[erad mg=%.2f L=%d] WALL reached at %d/%d; checkpointed, exit for resume.\n", mg, L, k, M); exit(0)
        end
    end
end

Jmap = Matrix(reduce(hcat, Jc)'); Hmap = Matrix(reduce(hcat, Hd)'); EFmap = Matrix(reduce(hcat, Ef)')
jldsave(outfile; mg=mg, L=L, ag=ag, N=N, c=c, T=T, dt=dt, dtheta=dtheta, links=[l1,l2],
        t=heat_t, Jmap=Jmap, Hmap=Hmap, EFmap=EFmap, hvac=hvac)
@printf("[erad mg=%.2f L=%d] DONE  injected E(t=0)=%.4f  -> %s\n", mg, L, sum(Hmap[1,:] .- hvac), outfile)

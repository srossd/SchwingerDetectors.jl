# Energy radiated after inserting a ±1/2 static probe-charge pair at separation L.
#
# θ=0 vacuum -> impose θ2π=+dtheta (=1/2) step over L links about the centre -> evolve under the
# stepped Hamiltonian. Radiated energy past a detector: E_rad = ∫₀ᵀ ⟨𝒥(x_D,t)⟩ dt.
#
# MEASUREMENT COST (confirmed by erad_meas_timing on the cluster stack): the FINITE-lattice
# energycurrents/energy_densities loop a global expectation per site -> O(N²) (~435/455 s per
# call at N=256!), while electricfields/chargecurrents use the fast local path (<0.2 s), and a
# SINGLE-bond EnergyCurrent expectation is O(N) (~2 s). So: detector 𝒥 via single-bond
# expectation, E-field/charge maps via the fast vectorized fns, NO full energycurrents/
# energy_densities. All measured every MEAS_STRIDE steps.
# (Cluster Schwinger lacks rehost()/standard_densities(): rehost by hand; vectorized fns only.)
#
#   julia --project=production production/erad_run.jl <m_over_g> <L_sites> <outdir>
using Statistics, Printf, JLD2, LinearAlgebra
using Schwinger, SchwingerDetectors
include(joinpath(@__DIR__, "persist.jl"))
mgtag(x) = replace(string(x), "." => "p")

mg = parse(Float64, ARGS[1]); L = parse(Int, ARGS[2]); outdir = get(ARGS, 3, "data/erad"); mkpath(outdir)
N = parse(Int, get(ENV, "N", "384")); ag = parse(Float64, get(ENV, "AG", "0.2"))
T = parse(Float64, get(ENV, "T", "50.0")); dt = parse(Float64, get(ENV, "DT", "0.2"))
maxbond = parse(Int, get(ENV, "MAXBOND", "128")); dtheta = parse(Float64, get(ENV, "DTHETA", "0.5"))
DETS = [parse(Int, x) for x in split(get(ENV, "DETS", "75,100"), ",")]   # detector offsets (sites); 75,100 = phys 15,20
MEAS_STRIDE = parse(Int, get(ENV, "MEAS_STRIDE", "5"))
PRINT_EVERY = parse(Int, get(ENV, "PRINT_EVERY", "5")); CKPT_EVERY = parse(Int, get(ENV, "CKPT_EVERY", "25"))
WALL = parse(Float64, get(ENV, "WALL_SECONDS", "61200"))
M = round(Int, T / dt); c = N ÷ 2
outfile = joinpath(outdir, "erad_mg$(mgtag(mg))_L$(L).jld2")
isfile(outfile) && (println("already done: $outfile"); exit(0))
ckd = joinpath(outdir, "erad_ckpt_mg$(mgtag(mg))_L$(L)"); mkpath(ckd)

m0 = build_model(N=N, F=1, q=1, ag=ag, mg=mg, theta2pi=0.0)
gsf = joinpath(outdir, "erad_gs_mg$(mgtag(mg)).bin")
gs = isfile(gsf) ? load_state(gsf, m0.H) : (g = groundstate(m0.H; energy_tol=1e-7, bonddim=96); save_state(gsf, g); g)

l1 = c - L ÷ 2; l2 = l1 + L
θ = fill(0.0, N); for n in l1:(l2 - 1); θ[n] += dtheta; end
latq = Lattice(N; F=1, q=1, a=ag, m=mg, θ2π=θ); Hq = Hamiltonian(latq; backend=:MPSKit)
detbonds = sort(unique([c + D for D in DETS if 1 <= c+D <= N-1] ∪ [c - D for D in DETS if 1 <= c-D <= N-1]))
jEops = [EnergyCurrent(latq, b; backend=:MPSKit) for b in detbonds]
detjE(s) = Float64[real(expectation(op, s)) for op in jEops]                 # O(N) per detector
efrow(s) = vec(real.(electricfields(s))); jcrow(s) = vec(real.(chargecurrents(s)))   # fast local

metaf = joinpath(ckd, "meta.jld2")
if isfile(metaf)
    d = load(metaf); k = d["k"]::Int
    t_meas = d["t_meas"]; DetJ = d["DetJ"]; EFmaps = d["EFmaps"]; JCmaps = d["JCmaps"]
    psi = load_state(joinpath(ckd, "psi.bin"), Hq)
    @printf("[erad mg=%.2f L=%d] RESUME at %d/%d\n", mg, L, k, M); flush(stdout)
else
    psi = MPSKitState(Hq, gs.psi)
    t_meas = [0.0]; DetJ = [detjE(psi)]; EFmaps = [efrow(psi)]; JCmaps = [jcrow(psi)]; k = 0
    @printf("[erad mg=%.2f L=%d] FRESH: θ=0 <E>=%.4f; string %d…%d (gx=%.2f) dθ=%.2f; N=%d T=%d M=%d maxbond=%d dets=%s meas_stride=%d\n",
            mg, L, mean(real.(electricfields(gs))), l1, l2, L*ag, dtheta, N, Int(T), M, maxbond, string(detbonds), MEAS_STRIDE); flush(stdout)
end

function checkpoint(k)
    save_state(joinpath(ckd, "psi.bin"), psi)
    tmp = metaf * ".tmp"; jldsave(tmp; k=k, t_meas=t_meas, DetJ=DetJ, EFmaps=EFmaps, JCmaps=JCmaps); mv(tmp, metaf; force=true)
end

t0 = time()
while k < M
    global k += 1
    global psi = evolve(psi, dt; nsteps=1, two_site=true, maxlinkdim=maxbond)[1]
    if k % MEAS_STRIDE == 0
        push!(t_meas, k*dt); push!(DetJ, detjE(psi)); push!(EFmaps, efrow(psi)); push!(JCmaps, jcrow(psi))
    end
    if k % PRINT_EVERY == 0
        el = time() - t0
        @printf("[erad mg=%.2f L=%d] step %d/%d  elapsed=%.2fh  proj-total=%.1fh\n", mg, L, k, M, el/3600, el/k*M/3600); flush(stdout)
    end
    if k % CKPT_EVERY == 0 || (time()-t0) > WALL
        checkpoint(k)
        (time()-t0) > WALL && k < M && (@printf("[erad mg=%.2f L=%d] WALL at %d/%d; checkpointed, exit.\n", mg, L, k, M); exit(0))
    end
end

DetJmat = Matrix(reduce(hcat, DetJ)')          # (nt_meas × ndet)
EFmap = Matrix(reduce(hcat, EFmaps)'); JCmap = Matrix(reduce(hcat, JCmaps)')
jldsave(outfile; mg=mg, L=L, ag=ag, N=N, c=c, T=T, dt=dt, dtheta=dtheta, links=[l1,l2], maxbond=maxbond,
        detbonds=detbonds, t=t_meas, DetJ=DetJmat, EFmap=EFmap, JCmap=JCmap)
@printf("[erad mg=%.2f L=%d] DONE maxbond=%d -> %s\n", mg, L, maxbond, outfile)

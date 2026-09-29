# Fast E-field / current heatmaps for one θ/π: evolve ONLY the physical state ψ (one state,
# ~28h at N=512) recording spacetime density maps. Resumable (checkpoint ψ + arrays), though it
# comfortably fits one 96h job. Much cheaper than the 5-state correlator, so heatmaps land in ~a day.
#   julia --project=production production/run_heat.jl <theta_over_pi> <outdir>
using Statistics, Printf, JLD2
using Schwinger, SchwingerDetectors
include(joinpath(@__DIR__, "config.jl"))
include(joinpath(@__DIR__, "persist.jl"))

theta_over_pi = parse(Float64, ARGS[1])
outdir = get(ARGS, 2, joinpath(@__DIR__, "..", "data", "qq_sweep"))
WALL = parse(Float64, get(ENV, "WALL_SECONDS", "334800"))
CKPT_EVERY = parse(Int, get(ENV, "CKPT_EVERY", "20"))

t2p = theta2pi(theta_over_pi); m = model(t2p); M = round(Int, T / DT)
tag = "th$(thetatag(theta_over_pi))"
ckd = joinpath(outdir, "heatckpt_$(tag)"); mkpath(ckd)
donef = joinpath(outdir, "heat_$(tag).jld2")
isfile(donef) && (println("[$tag heat] already DONE"); exit(0))
mb = MAXBOND
fwd(s) = evolve(s, DT; nsteps=SUBSTEP, two_site=true, maxlinkdim=mb)[1]
erow(s) = vec(real.(electricfields(s))); jrow(s) = vec(real.(chargecurrents(s)))

metaf = joinpath(ckd, "meta.jld2")
if isfile(metaf)
    d = load(metaf); k = d["k"]::Int
    heat_t = d["heat_t"]; Emaps = d["Emaps"]; Jmaps = d["Jmaps"]
    psi = load_state(joinpath(ckd,"psi.bin"), m.H)
    @printf("[%s heat] RESUME at %d/%d\n", tag, k, M); flush(stdout)
else
    gs = load_state(joinpath(outdir, "gs_theta$(thetatag(theta_over_pi)).bin"), m.H)
    psi = wilson_line_quench(m, wilson_endpoints()...; gs=gs).state
    heat_t = [0.0]; Emaps = [erow(psi)]; Jmaps = [jrow(psi)]; k = 0
    @printf("[%s heat] FRESH: N=%d ag=%.2f θ2π=%.3f T=%.0f M=%d dt=%.2f\n", tag, N, AG, t2p, T, M, DT); flush(stdout)
end

t0 = time()
while k < M
    global k += 1
    global psi = fwd(psi)
    if k % HEAT_STRIDE == 0; push!(heat_t, k*DT); push!(Emaps, erow(psi)); push!(Jmaps, jrow(psi)); end
    if k % CKPT_EVERY == 0 || (time()-t0) > WALL
        save_state(joinpath(ckd,"psi.bin"), psi)
        tmp = metaf*".tmp"; jldsave(tmp; k=k, heat_t=heat_t, Emaps=Emaps, Jmaps=Jmaps); mv(tmp, metaf; force=true)
        @printf("[%s heat] step %d/%d elapsed=%.1fh\n", tag, k, M, (time()-t0)/3600); flush(stdout)
        (time()-t0) > WALL && k < M && (println("[$tag heat] WALL reached; checkpointed, exit for resume."); exit(0))
    end
end

Emap = Matrix(reduce(hcat, Emaps)'); Jmap = Matrix(reduce(hcat, Jmaps)')
jldsave(donef; theta_over_pi=theta_over_pi, t2p=t2p, N=N, ag=AG, mg=MG, center=CENTER,
        wilson=collect(wilson_endpoints()), Rlist=RLIST, T=T, dt=DT, maxbond=mb,
        heat_t=heat_t, Emap=Emap, Jmap=Jmap)
@printf("[%s heat] DONE -> %s\n", tag, donef)

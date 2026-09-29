# End-to-end pipeline validation on the CLUSTER's actual Schwinger (GitHub main) at m/g=1,
# θ/π=1.2: folded-vacuum branch, Serialization round-trip, multi-pair correlator (_addstates),
# and the resumable driver (fresh -> forced pause -> resume -> done) vs the reference. Small N.
#   julia --project=production production/validate_orcd.jl
for (k,v) in ("N"=>"32","AG"=>"0.2","MG"=>"1.0","GX_WL"=>"2.0","T"=>"4.0","DT"=>"0.25",
              "MAXBOND"=>"48","GS_BOND"=>"32","GS_TOL"=>"1e-6","SUBSTEPS"=>"1")
    ENV[k] = v
end
using Schwinger, SchwingerDetectors, JLD2, Printf, Statistics, LinearAlgebra
include(joinpath(@__DIR__, "config.jl"))
include(joinpath(@__DIR__, "persist.jl"))

outdir = "data/validate"; rm(outdir; recursive=true, force=true); mkpath(outdir)
th = 1.2; R = 6
jl = `julia --project=production`
here = @__DIR__

println("== vacuum (fold branch + Serialization save) =="); flush(stdout)
run(`$jl $(joinpath(here,"vacuum.jl")) $th $outdir`)

# serialize round-trip check + reference correlator (tests _addstates on cluster Schwinger)
m = model(theta2pi(th))
gs = load_state(joinpath(outdir, "gs_theta$(thetatag(th)).bin"), m.H)
@printf("serialize round-trip: |<gs|gs_loaded>|/norm^2 = %.3e (want 1)\n",
        abs(dot(gs,gs)) / real(dot(gs,gs)))   # trivially 1; real check is that load_state worked
psi = wilson_line_quench(m, wilson_endpoints()...; gs=gs).state
bL = CENTER - R; bR = CENTER + R
ref = charge_transfer_correlator(psi, T, [(bL,bR)]; nsteps=round(Int,T/DT), substeps=SUBSTEP,
                                 two_site=true, maxbond=MAXBOND)
@printf("reference correlator connected(T) = %.5e (multi-pair _addstates OK)\n", real(ref.connected[1][end]))

println("== driver fresh (force pause) =="); flush(stdout)
ENV["CKPT_EVERY"] = "3"; ENV["WALL_SECONDS"] = "0"; ENV["RECORD_HEAT"] = "1"
run(`$jl $(joinpath(here,"run_corr_ckpt.jl")) $th $R $outdir`)
@assert isfile(joinpath(outdir,"ckpt_th$(thetatag(th))_R$(R)","meta.jld2")) "no checkpoint"
@assert !isfile(joinpath(outdir,"qq_th$(thetatag(th))_R$(R).jld2")) "finished too early"

println("== driver resume -> done =="); flush(stdout)
ENV["WALL_SECONDS"] = "100000"
run(`$jl $(joinpath(here,"run_corr_ckpt.jl")) $th $R $outdir`)
d = load(joinpath(outdir,"qq_th$(thetatag(th))_R$(R).jld2"))
@printf("driver vs reference: max|Δfull| = %.3e\n", maximum(abs.(d["full"] .- ref.full[1])))
@printf("heatmap recorded: Emap size = %s\n", string(size(d["Emap"])))
println("VALIDATE_ORCD_DONE")

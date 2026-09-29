# Aggregation + plots for the per-(θ,R) production outputs (qq_th<tag>_R<R>.jld2).
# Produces:
#   (1) <QQ> vs θ/π for each R      -> qq_vs_theta.png (+ full-magnitude variant)
#   (2) E-field & current heatmaps  -> heat_E_theta<tag>.png, heat_J_theta<tag>.png (per θ)
#       (heatmaps come from the RECORD_HEAT task, i.e. the smallest-R file per θ)
#
#   julia --project=production production/plots.jl <outdir>
ENV["GKSwstype"] = "100"   # headless GR
using JLD2, Printf, Plots
gr()

outdir = get(ARGS, 1, joinpath(@__DIR__, "..", "data", "qq_sweep"))
plotdir = joinpath(outdir, "plots"); mkpath(plotdir)

files = filter(f -> occursin(r"^qq_th.*_R\d+\.jld2$", basename(f)), readdir(outdir; join=true))
isempty(files) && error("no qq_th*_R*.jld2 in $outdir")
recs = [load(f) for f in files]

thetas = sort(unique(r["theta_over_pi"] for r in recs))
Rs     = sort(unique(r["R"] for r in recs))
ag     = recs[1]["ag"]
getrec(th, R) = recs[findfirst(r -> r["theta_over_pi"]==th && r["R"]==R, recs)]
@printf("θ/π = %s ;  R = %s (phys %s)\n", string(thetas), string(Rs), string(Rs .* ag))

# ---- Plot 1: <QQ> vs θ/π, one series per R ----
plt = plot(xlabel="θ/π", ylabel="Re ⟨Q(-R) Q(R)⟩_c", title="Charge-transfer correlator vs θ  (m/g=1, ag=$ag)", legend=:best)
for R in Rs
    ys = [real(getrec(th, R)["connected"][end]) for th in thetas]
    plot!(plt, thetas, ys; label="R=$R (x=$(R*ag))", marker=:circle)
end
vline!(plt, [1.0]; ls=:dash, color=:gray, label="θ=π")
savefig(plt, joinpath(plotdir, "qq_vs_theta.png"))

pltf = plot(xlabel="θ/π", ylabel="|⟨Q(-R) Q(R)⟩|", title="Charge-transfer correlator (full) vs θ", legend=:best)
for R in Rs
    plot!(pltf, thetas, [abs(getrec(th, R)["full"][end]) for th in thetas]; label="R=$R", marker=:circle)
end
vline!(pltf, [1.0]; ls=:dash, color=:gray, label="θ=π")
savefig(pltf, joinpath(plotdir, "qq_full_vs_theta.png"))
println("wrote qq_vs_theta.png, qq_full_vs_theta.png")

# ---- Plot 2: E-field & current heatmaps per θ (prefer dedicated heat_th*.jld2 files) ----
for th in thetas
    tagf = replace(string(th), "." => "p")
    heatfile = joinpath(outdir, "heat_th$(tagf).jld2")
    r = nothing
    if isfile(heatfile)
        r = load(heatfile)
    else
        for R in Rs
            rr = getrec(th, R)
            if haskey(rr, "Emap") && !isempty(rr["Emap"]); r = rr; break; end
        end
    end
    r === nothing && (@warn "no heatmap data for θ/π=$th"; continue)
    tag = replace(string(th), "." => "p")
    heat_t = r["heat_t"]; Emap = r["Emap"]; Jmap = r["Jmap"]
    nbond = size(Emap, 2)
    x = ((1:nbond) .- r["center"]) .* ag
    for (Z, nm, ttl) in ((Emap, "E", "Electric field"), (Jmap, "J", "Charge current j¹"))
        cmax = maximum(abs, Z); cmax = cmax == 0 ? 1.0 : cmax
        h = heatmap(x, heat_t, Z; xlabel="x (g·x, rel. centre)", ylabel="t (g·t)",
                    title="$ttl,  θ/π=$th", color=:balance, clims=(-cmax, cmax))
        for R in Rs; vline!(h, [-R*ag, R*ag]; ls=:dot, color=:black, label="", alpha=0.4); end
        savefig(h, joinpath(plotdir, "heat_$(nm)_theta$(tag).png"))
    end
end
println("wrote heatmaps; plots in $plotdir")

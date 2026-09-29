# Aggregation + plots. Reads all qq_theta*.jld2 in <outdir> and produces:
#   (1) <QQ> vs θ/π for each R      -> qq_vs_theta.png
#   (2) E-field & current heatmaps  -> heat_E_theta<..>.png, heat_J_theta<..>.png  (per θ)
#
#   julia --project=production production/plots.jl <outdir>
ENV["GKSwstype"] = "100"   # headless GR
using JLD2, Printf, Plots
gr()

outdir = get(ARGS, 1, joinpath(@__DIR__, "..", "data", "qq_sweep"))
plotdir = joinpath(outdir, "plots"); mkpath(plotdir)

files = sort(filter(f -> startswith(basename(f), "qq_theta") && endswith(f, ".jld2"),
                    readdir(outdir; join = true)))
isempty(files) && error("no qq_theta*.jld2 found in $outdir")

data = [load(f) for f in files]
order = sortperm([d["theta_over_pi"] for d in data])
data = data[order]
thetas = [d["theta_over_pi"] for d in data]
Rlist = data[1]["Rlist"]
ag = data[1]["ag"]
@printf("loaded %d θ points: %s;  R = %s\n", length(data), string(thetas), string(Rlist))

# ---- Plot 1: <QQ> vs θ/π for each R (connected, real part) ----
plt = plot(xlabel = "θ/π", ylabel = "Re ⟨Q(-R) Q(R)⟩_c", title = "Charge-transfer correlator vs θ",
           legend = :best, marker = :circle)
for (ri, R) in enumerate(Rlist)
    ys = [real(d["connected"][ri][end]) for d in data]
    plot!(plt, thetas, ys; label = "R=$R", marker = :circle)
end
vline!(plt, [1.0]; ls = :dash, color = :gray, label = "θ=π (transition)")
savefig(plt, joinpath(plotdir, "qq_vs_theta.png"))
# also the full (unsubtracted) correlator magnitude, for reference
plt2 = plot(xlabel = "θ/π", ylabel = "|⟨Q(-R) Q(R)⟩|", title = "Charge-transfer correlator (full) vs θ",
            legend = :best)
for (ri, R) in enumerate(Rlist)
    plot!(plt2, thetas, [abs(d["full"][ri][end]) for d in data]; label = "R=$R", marker = :circle)
end
vline!(plt2, [1.0]; ls = :dash, color = :gray, label = "θ=π")
savefig(plt2, joinpath(plotdir, "qq_full_vs_theta.png"))
println("wrote qq_vs_theta.png, qq_full_vs_theta.png")

# ---- Plot 2: E-field and current heatmaps, per θ ----
for d in data
    tag = replace(string(d["theta_over_pi"]), "." => "p")
    heat_t = d["heat_t"]; Emap = d["Emap"]; Jmap = d["Jmap"]
    nbond = size(Emap, 2)
    x = ((1:nbond) .- d["center"]) .* ag        # physical position relative to centre
    for (Z, name, cl, ttl) in ((Emap, "E", :balance, "Electric field"),
                               (Jmap, "J", :balance, "Charge current j¹"))
        cmax = maximum(abs, Z); cmax = cmax == 0 ? 1.0 : cmax
        h = heatmap(x, heat_t, Z; xlabel = "x (relative to centre, g·x)", ylabel = "t (g·t)",
                    title = "$ttl,  θ/π=$(d["theta_over_pi"])", color = cl, clims = (-cmax, cmax))
        # mark detector positions
        for R in Rlist; vline!(h, [-R*ag, R*ag]; ls = :dot, color = :black, label = "", alpha = 0.4); end
        savefig(h, joinpath(plotdir, "heat_$(name)_theta$(tag).png"))
    end
end
println("wrote heatmaps for θ/π = ", string(thetas))
println("plots in $plotdir")

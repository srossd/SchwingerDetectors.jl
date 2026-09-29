# Aggregate the energy-radiation runs (erad_mg<>_L<>.jld2):
#   (1) radiated energy E_rad(L) = ∫₀ᵀ ⟨𝒥(x_D,t)⟩ dt at a fixed detector, per m/g
#   (2) E-field and charge-current spacetime heatmaps per (m/g, L)
# Tests: E_rad should saturate in L at m/g=0 (screening) but rise ~linearly at m/g=1 (confinement).
#   julia --project=production production/erad_plots.jl <outdir>
ENV["GKSwstype"] = "100"
using JLD2, Printf, Plots
gr()

outdir = get(ARGS, 1, "data/erad"); plotdir = joinpath(outdir, "plots"); mkpath(plotdir)
files = filter(f -> occursin(r"^erad_mg.*_L\d+\.jld2$", basename(f)), readdir(outdir; join=true))
isempty(files) && error("no erad_mg*_L*.jld2 in $outdir")
recs = [load(f) for f in files]
masses = sort(unique(r["mg"] for r in recs)); Ls = sort(unique(r["L"] for r in recs))
ag = recs[1]["ag"]; c = recs[1]["c"]
getrec(mg, L) = recs[findfirst(r -> r["mg"]==mg && r["L"]==L, recs)]
trapz(t, y) = sum((@view(y[2:end]) .+ @view(y[1:end-1]))/2 .* diff(t))
@printf("masses=%s  Ls=%s (phys %s)\n", string(masses), string(Ls), string(Ls .* ag))
DET = parse(Int, get(ENV, "DET_OFF", "100"))   # detector offset (sites); 100 = phys 20

# ---- (1) radiated energy vs L ----
plt = plot(xlabel="L (sites)", ylabel="E_rad = ∫𝒥(x_D,t) dt",
           title="Radiated energy vs L  (±1/2 charges, x_D=+$(round(DET*ag;digits=1)))", legend=:topleft)
for mg in masses
    ys = Float64[]
    for L in Ls
        r = getrec(mg, L); col = findfirst(==(c+DET), r["detbonds"])
        push!(ys, col === nothing ? NaN : trapz(r["t"], r["DetJ"][:, col]))
    end
    plot!(plt, Ls, ys; label="m/g=$mg", marker=:circle)
    @printf("E_rad(m/g=%.1f): %s\n", mg, string(round.(ys; sigdigits=4)))
end
savefig(plt, joinpath(plotdir, "erad_vs_L.png"))

# ---- (2) E-field & charge-current heatmaps per (m/g, L) ----
for r in recs
    for (Z, nm, ttl) in ((r["EFmap"], "E", "Electric field"), (r["JCmap"], "jq", "Charge current j¹"))
        nb = size(Z, 2); x = ((1:nb) .- c) .* ag; cmax = maximum(abs, Z); cmax = cmax==0 ? 1.0 : cmax
        h = heatmap(x, r["t"], Z; xlabel="x (g·x)", ylabel="t (g·t)", color=:balance, clims=(-cmax, cmax),
                    title="$ttl  m/g=$(r["mg"]) L=$(r["L"]) (gx=$(round(r["L"]*ag;digits=1)))")
        savefig(h, joinpath(plotdir, "heat_$(nm)_mg$(replace(string(r["mg"]),"."=>"p"))_L$(r["L"]).png"))
    end
end
println("plots in $plotdir")

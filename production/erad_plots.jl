# Aggregate the energy-radiation runs (erad_mg<>_L<>.jld2):
#   (1) radiated energy E_rad(L) = ∫₀ᵀ ⟨𝒥(x_D,t)⟩ dt at a fixed detector, per m/g
#   (2) injected energy Σ(h−h_vac) at t=0, per m/g
#   (3) energy-current heatmaps ⟨𝒥(x,t)⟩ per (m/g, L)
# Tests: E should saturate in L at m/g=0 (screening) but rise ~linearly at m/g=1 (confinement).
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

DET = parse(Int, get(ENV, "DET_OFF", "25"))   # detector at bond c+DET (physical DET*ag), outside all strings

# ---- (1) radiated energy vs L ----
plt = plot(xlabel="L (sites)", ylabel="E_rad = ∫𝒥(x_D,t) dt",
           title="Radiated energy vs L  (±1/2 charges, x_D=+$(round(DET*ag;digits=1)))", legend=:topleft)
erad = Dict()
for mg in masses
    ys = Float64[]
    for L in Ls
        r = getrec(mg, L); t = r["t"]; Jmap = r["Jmap"]
        col = min(c + DET, size(Jmap, 2))
        push!(ys, trapz(t, Jmap[:, col]))
    end
    erad[mg] = ys
    plot!(plt, Ls, ys; label="m/g=$mg", marker=:circle)
    @printf("E_rad(m/g=%.1f): %s\n", mg, string(round.(ys; sigdigits=4)))
end
savefig(plt, joinpath(plotdir, "erad_vs_L.png"))

# ---- (2) injected energy vs L ----
plti = plot(xlabel="L (sites)", ylabel="Σ(h−h_vac) at t=0  (injected energy)",
            title="Injected energy vs L (±1/2 charges)", legend=:topleft)
for mg in masses
    ys = [sum(getrec(mg, L)["Hmap"][1, :] .- getrec(mg, L)["hvac"]) for L in Ls]
    plot!(plti, Ls, ys; label="m/g=$mg", marker=:square)
    @printf("E_inject(m/g=%.1f): %s\n", mg, string(round.(ys; sigdigits=4)))
end
savefig(plti, joinpath(plotdir, "einject_vs_L.png"))

# ---- (3) energy-current heatmaps ----
for r in recs
    Jmap = r["Jmap"]; t = r["t"]; nb = size(Jmap, 2)
    x = ((1:nb) .- c) .* ag
    cmax = maximum(abs, Jmap); cmax = cmax == 0 ? 1.0 : cmax
    h = heatmap(x, t, Jmap; xlabel="x (g·x)", ylabel="t (g·t)", color=:balance, clims=(-cmax, cmax),
                title="⟨𝒥(x,t)⟩  m/g=$(r["mg"]) L=$(r["L"]) (gx=$(round(r["L"]*ag;digits=1)))")
    savefig(h, joinpath(plotdir, "heat_J_mg$(replace(string(r["mg"]),"."=>"p"))_L$(r["L"]).png"))
end
println("plots in $plotdir")

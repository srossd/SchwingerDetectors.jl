# QQ (Wilson-line quench) E-field & charge-current heatmaps, up to whatever time the run_heat
# jobs have reached — reads the finished heat_th<tag>.jld2 if present, else the (atomically
# written) checkpoint heatckpt_th<tag>/meta.jld2. Applies 1-2-1 spatial smoothing to remove the
# staggered even/odd sawtooth.
#   julia --project=production production/qq_heat_plots.jl [outdir]
ENV["GKSwstype"] = "100"
using JLD2, Printf, Plots
gr()

outdir = get(ARGS, 1, "data/qq_sweep")
N = parse(Int, get(ENV, "N", "512")); ag = parse(Float64, get(ENV, "AG", "0.2")); c = N ÷ 2
plotdir = joinpath(outdir, "plots"); mkpath(plotdir)

# one pass of 1-2-1 (binomial) smoothing along the spatial axis (columns) of an nt×nx map
function smooth121(Z)
    S = copy(Z)
    @views for j in 2:size(Z, 2)-1
        S[:, j] .= (Z[:, j-1] .+ 2 .* Z[:, j] .+ Z[:, j+1]) ./ 4
    end
    return S
end

for th in (0.6, 0.8, 1.0, 1.2, 1.4)
    tag = replace(string(th), "." => "p")
    donef = joinpath(outdir, "heat_th$(tag).jld2")
    ckf = joinpath(outdir, "heatckpt_th$(tag)", "meta.jld2")
    heat_t = nothing; Emap = nothing; Jmap = nothing; src = ""
    if isfile(donef)
        d = load(donef); heat_t = d["heat_t"]; Emap = d["Emap"]; Jmap = d["Jmap"]; src = "done"
    elseif isfile(ckf)
        d = load(ckf); heat_t = d["heat_t"]
        Emap = Matrix(reduce(hcat, d["Emaps"])'); Jmap = Matrix(reduce(hcat, d["Jmaps"])'); src = "checkpoint(k=$(d["k"]))"
    else
        @warn "no heat data for θ/π=$th"; continue
    end
    tmax = heat_t[end]
    for (Z0, nm, ttl) in ((Emap, "E", "Electric field"), (Jmap, "jq", "Charge current j¹"))
        Z = smooth121(Z0); nb = size(Z, 2); x = ((1:nb) .- c) .* ag
        cmax = maximum(abs, Z); cmax = cmax == 0 ? 1.0 : cmax
        h = heatmap(x, heat_t, Z; xlabel="x (g·x, rel. centre)", ylabel="t (g·t)", color=:balance,
                    clims=(-cmax, cmax), title="$ttl (1-2-1 smoothed)  θ/π=$th  (t ≤ $(round(tmax;digits=1)))")
        savefig(h, joinpath(plotdir, "qqheat_$(nm)_theta$(tag).png"))
    end
    @printf("θ/π=%.1f [%s]: t≤%.1f, %d frames, %d bonds -> qqheat_{E,jq}_theta%s.png\n",
            th, src, tmax, length(heat_t), size(Emap, 2), tag)
end
println("QQ heatmaps written to $plotdir")

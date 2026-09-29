# Cumulative radiated energy ∫₀ᵀ 𝒥(x_D,t) dt vs upper limit T, for both m/g at a fixed L.
# Detector offsets (sites from centre) via ENV DETS (default 75,100 = physical 15,20).
#   julia --project=production production/erad_cumulative.jl <outdir> <L_sites> [out.png]
ENV["GKSwstype"] = "100"
using JLD2, Printf, Plots
gr()

outdir = ARGS[1]; L = parse(Int, ARGS[2])
outpng = get(ARGS, 3, joinpath(outdir, "plots", "cumulative_L$(L).png")); mkpath(dirname(outpng))
dets = [parse(Int, x) for x in split(get(ENV, "DETS", "75,100"), ",")]

plt = plot(xlabel="T (g·t)", ylabel="∫₀ᵀ ⟨𝒥(x_D,t)⟩ dt", legend=:topleft,
           title="Cumulative radiated energy vs T  (±1/2 charges, L=$L sites)")
found = false
for mg in (0.0, 1.0)
    f = joinpath(outdir, "erad_mg$(replace(string(mg),"."=>"p"))_L$(L).jld2")
    isfile(f) || (@warn "missing $f"; continue)
    global found = true
    d = load(f); t = d["t"]; Jmap = d["Jmap"]; c = d["c"]; ag = d["ag"]
    ls = mg == 0.0 ? :solid : :dash
    for D in dets
        col = min(c + D, size(Jmap, 2)); j = Jmap[:, col]
        C = zeros(length(t)); for k in 2:length(t); C[k] = C[k-1] + (j[k]+j[k-1])/2*(t[k]-t[k-1]); end
        plot!(plt, t, C; label="m/g=$mg, x_D=$(round(D*ag; digits=1))", ls=ls, lw=2)
        @printf("m/g=%.1f x_D=%.1f: E_rad(T=%.0f)=%.4e\n", mg, D*ag, t[end], C[end])
    end
end
found || error("no erad files for L=$L in $outdir")
savefig(plt, outpng); println("wrote $outpng")

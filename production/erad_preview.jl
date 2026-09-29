# Quick 3-panel preview of one energy-radiation run: electric field, excess energy density,
# and energy current 𝒥, as spacetime heatmaps. (JLD2 + Plots only — no Schwinger, loads fast.)
#   julia --project=production production/erad_preview.jl <erad_*.jld2> <out.png>
ENV["GKSwstype"] = "100"
using JLD2, Printf, Plots
gr()

f = ARGS[1]; outpng = ARGS[2]
d = load(f); c = d["c"]; ag = d["ag"]; t = d["t"]
Jmap = d["Jmap"]; Hmap = d["Hmap"]; EFmap = d["EFmap"]; hvac = d["hvac"]
hexc = Hmap .- reshape(hvac, 1, :)                    # energy density above θ=0 vacuum
xof(Z) = ((1:size(Z, 2)) .- c) .* ag
sym(Z) = (m = maximum(abs, Z); m == 0 ? 1.0 : m)
pan(Z, ttl) = heatmap(xof(Z), t, Z; title=ttl, xlabel="x (g·x)", ylabel="t (g·t)",
                      color=:balance, clims=(-sym(Z), sym(Z)))
plt = plot(pan(EFmap, "Electric field"), pan(hexc, "Energy density (excess)"),
           pan(Jmap, "Energy current 𝒥"); layout=(1, 3), size=(1350, 400),
           plot_title="±1/2 charges  m/g=$(d["mg"]) L=$(d["L"]) (N=$(d["N"]), gx=$(round(d["L"]*ag; digits=2)))")
savefig(plt, outpng)
@printf("injected E (t=0)=%.4f;  wrote %s\n", sum(hexc[1, :]), outpng)

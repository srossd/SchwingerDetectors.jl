# Stage 2 (expensive): for one θ/π, load the verified vacuum, apply the Wilson-line quench,
# evolve to T, and in a SINGLE streaming pass compute (a) the charge-transfer correlators
# ⟨Q(c-R)Q(c+R)⟩ for all R, and (b) the electric-field / current spacetime maps (for heatmaps).
# A cheap low-resolution method=:both cross-check is run first as a correctness gate.
#
#   julia --project=production production/run_theta.jl <theta_over_pi> <outdir>
using Statistics, Printf, JLD2, LinearAlgebra
include(joinpath(@__DIR__, "config.jl"))

theta_over_pi = parse(Float64, ARGS[1])
outdir = get(ARGS, 2, joinpath(@__DIR__, "..", "data", "qq_sweep"))
mkpath(outdir)
t2p = theta2pi(theta_over_pi)
m = model(t2p)
M = round(Int, T / DT)
@printf("[run] θ/π=%.3f θ2π=%.3f N=%d ag=%.3f m/g=%.3f  T=%.1f M=%d dt=%.3f maxbond=%d\n",
        theta_over_pi, t2p, N, AG, MG, T, M, DT, MAXBOND); flush(stdout)

# vacuum (from stage 1) + Wilson-line quench
gspath = joinpath(outdir, "gs_theta$(thetatag(theta_over_pi)).jld2")
gs = loadstate(gspath, m.H)
Eavg = mean(real.(electricfields(gs)))
@printf("[run] loaded vacuum <E_avg>=%.4f from %s\n", Eavg, gspath); flush(stdout)
l1, l2 = wilson_endpoints()
psi = wilson_line_quench(m, l1, l2; gs = gs).state
@printf("[run] wilson line sites %d…%d (gx=%.1f)\n", l1, l2, GX_WL); flush(stdout)

# --- correctness gate: low-res both-method cross-check on the closest detector ---
chk = charge_transfer_correlator(psi, min(T, 3.0), (CENTER - RLIST[1], CENTER + RLIST[1]);
                                 nsteps = 12, substeps = 2, two_site = true,
                                 maxbond = MAXBOND, method = :both)
@printf("[run] cross-check (R=%d, T≤3): max|accumulated-matrix| = %.3e\n", RLIST[1], chk.discrepancy)
flush(stdout)

# --- heatmap recorder: electric field & current every HEAT_STRIDE steps ---
heat_t = Float64[]; Emaps = Vector{Vector{Float64}}(); Jmaps = Vector{Vector{Float64}}()
function recorder(ψ, t, k)
    (k % HEAT_STRIDE == 0) || return
    push!(heat_t, t)
    push!(Emaps, real.(electricfields(ψ)))
    push!(Jmaps, real.(chargecurrents(ψ)))
    return
end

# --- production: streaming correlator for all R in one evolution ---
pairs = [(CENTER - R, CENTER + R) for R in RLIST]
t_run = @elapsed res = charge_transfer_correlator(psi, T, pairs; nsteps = M, substeps = SUBSTEP,
                                                  two_site = true, maxbond = MAXBOND,
                                                  observe = recorder)
@printf("[run] evolution done in %.1f min\n", t_run/60); flush(stdout)
for (i, R) in enumerate(RLIST)
    @printf("[run] R=%2d:  full(T)=(% .4e,% .4e)  connected(T)=% .4e\n",
            R, real(res.full[i][end]), imag(res.full[i][end]), real(res.connected[i][end]))
end

Emap = reduce(hcat, Emaps)' |> Matrix   # (n_heat × nbond)
Jmap = reduce(hcat, Jmaps)' |> Matrix
outfile = joinpath(outdir, "qq_theta$(thetatag(theta_over_pi)).jld2")
jldsave(outfile;
        theta_over_pi = theta_over_pi, t2p = t2p, N = N, ag = AG, mg = MG,
        center = CENTER, wilson = collect((l1, l2)), Rlist = RLIST, T = T, dt = DT,
        maxbond = MAXBOND, substeps = SUBSTEP, Eavg = Eavg,
        t = res.t, full = res.full, connected = res.connected, QL = res.QL, QR = res.QR,
        heat_t = heat_t, Emap = Emap, Jmap = Jmap,
        cross_check = chk.discrepancy)
println("[run] wrote $outfile")

# Resumable charge-transfer correlator for ONE (θ/π, R) detector pair.
#
# Implements the validated streaming :accumulated recursion (identical to
# SchwingerDetectors.charge_transfer_correlator) with checkpoint/resume so a T=90 evolution
# can span several ≤96h jobs. Optionally records E-field/current spacetime maps (RECORD_HEAT=1)
# for the heatmaps in the same pass. On wall-time budget exhaustion it checkpoints and exits 0;
# the next (dependent) job resumes from the checkpoint. When it reaches M it writes the final
# result + a DONE marker and later jobs no-op.
#
#   julia --project=production production/run_corr_ckpt.jl <theta_over_pi> <R> <outdir>
using Statistics, Printf, JLD2, LinearAlgebra
using Schwinger, SchwingerDetectors
const CTC = SchwingerDetectors
include(joinpath(@__DIR__, "config.jl"))

theta_over_pi = parse(Float64, ARGS[1])
R      = parse(Int, ARGS[2])
outdir = get(ARGS, 3, joinpath(@__DIR__, "..", "data", "qq_sweep"))
WALL   = parse(Float64, get(ENV, "WALL_SECONDS", "334800"))  # ~93h, margin under the 96h limit
CKPT_EVERY = parse(Int, get(ENV, "CKPT_EVERY", "15"))
RECORD_HEAT = get(ENV, "RECORD_HEAT", "0") == "1"

t2p = theta2pi(theta_over_pi)
m   = model(t2p)
M   = round(Int, T / DT)
bL  = CENTER - R; bR = CENTER + R
tag = "th$(thetatag(theta_over_pi))_R$(R)"
ckd = joinpath(outdir, "ckpt_$(tag)"); mkpath(ckd)
donef = joinpath(outdir, "qq_$(tag).jld2")
isfile(donef) && (println("[$tag] already DONE ($donef); nothing to do."); exit(0))

jL = ChargeCurrent(m.lat, bL; backend=:MPSKit)
jR = ChargeCurrent(m.lat, bR; backend=:MPSKit)
mb = MAXBOND
fwd(s) = evolve(s, DT; nsteps=SUBSTEP, two_site=true, maxlinkdim=mb)[1]
jexp(j,s) = real(expectation(j, s))

# atomic savestate (write to tmp, rename)
function asave(path, st); tmp = path * ".tmp"; savestate(tmp, st); mv(tmp, path; force=true); end

metaf = joinpath(ckd, "meta.jld2")
if isfile(metaf)   # ---- resume ----
    d = load(metaf)
    k = d["k"]::Int
    full = d["full"]::Vector{ComplexF64}
    jLt = d["jLt"]::Vector{Float64}; jRt = d["jRt"]::Vector{Float64}
    heat_t = d["heat_t"]::Vector{Float64}
    Emaps = d["Emaps"]::Vector{Vector{Float64}}; Jmaps = d["Jmaps"]::Vector{Vector{Float64}}
    psi = loadstate(joinpath(ckd,"psi.jld2"), m.H)
    Rv  = loadstate(joinpath(ckd,"Rv.jld2"),  m.H); Ru  = loadstate(joinpath(ckd,"Ru.jld2"),  m.H)
    e0v = loadstate(joinpath(ckd,"e0v.jld2"), m.H); e0u = loadstate(joinpath(ckd,"e0u.jld2"), m.H)
    @printf("[%s] RESUME at step %d/%d\n", tag, k, M); flush(stdout)
else               # ---- fresh start ----
    gs = loadstate(joinpath(outdir, "gs_theta$(thetatag(theta_over_pi)).jld2"), m.H)
    psi = wilson_line_quench(m, wilson_endpoints()...; gs=gs).state
    a0 = act(jR, psi); b0 = act(jL, psi)
    Rv = DT*a0; Ru = DT*b0; e0v = a0; e0u = b0
    full = zeros(ComplexF64, M+1); jLt = zeros(M+1); jRt = zeros(M+1)
    jLt[1] = jexp(jL, psi); jRt[1] = jexp(jR, psi)
    heat_t = Float64[]; Emaps = Vector{Vector{Float64}}(); Jmaps = Vector{Vector{Float64}}()
    if RECORD_HEAT; push!(heat_t, 0.0); push!(Emaps, real.(electricfields(psi))); push!(Jmaps, real.(chargecurrents(psi))); end
    k = 0
    @printf("[%s] FRESH start: N=%d ag=%.2f m/g=%.1f θ2π=%.3f bonds=(%d,%d) T=%.0f M=%d dt=%.2f maxbond=%d\n",
            tag, N, AG, MG, t2p, bL, bR, T, M, DT, mb); flush(stdout)
end

function checkpoint(k)
    asave(joinpath(ckd,"psi.jld2"), psi); asave(joinpath(ckd,"Rv.jld2"), Rv)
    asave(joinpath(ckd,"Ru.jld2"), Ru);  asave(joinpath(ckd,"e0v.jld2"), e0v); asave(joinpath(ckd,"e0u.jld2"), e0u)
    tmp = metaf * ".tmp"
    jldsave(tmp; k=k, full=full, jLt=jLt, jRt=jRt, heat_t=heat_t, Emaps=Emaps, Jmaps=Jmaps,
            theta_over_pi=theta_over_pi, t2p=t2p, N=N, ag=AG, mg=MG, R=R, T=T, dt=DT, maxbond=mb, M=M)
    mv(tmp, metaf; force=true)
end

t0 = time()
while k < M
    global k += 1
    global psi = fwd(psi)
    ak = act(jR, psi); bk = act(jL, psi)
    global Rv = CTC._addstates(fwd(Rv), DT*ak); global e0v = fwd(e0v)
    global Ru = CTC._addstates(fwd(Ru), DT*bk); global e0u = fwd(e0u)
    vk = CTC._addstates(Rv, (-DT/2)*ak, (-DT/2)*e0v)
    uk = CTC._addstates(Ru, (-DT/2)*bk, (-DT/2)*e0u)
    full[k+1] = dot(uk, vk)
    jLt[k+1] = jexp(jL, psi); jRt[k+1] = jexp(jR, psi)
    if RECORD_HEAT && (k % HEAT_STRIDE == 0)
        push!(heat_t, k*DT); push!(Emaps, real.(electricfields(psi))); push!(Jmaps, real.(chargecurrents(psi)))
    end
    if k % CKPT_EVERY == 0 || (time()-t0) > WALL
        checkpoint(k)
        @printf("[%s] step %d/%d  full=%.4e%+.4eim  elapsed=%.1fh\n",
                tag, k, M, real(full[k+1]), imag(full[k+1]), (time()-t0)/3600); flush(stdout)
    end
    if (time()-t0) > WALL && k < M
        @printf("[%s] WALL budget reached at step %d/%d; checkpointed, exiting for resume.\n", tag, k, M)
        exit(0)
    end
end

# ---- finished: finalize + write result, clear checkpoint ----
QL = CTC._cumintegral(jLt, M, DT); QR = CTC._cumintegral(jRt, M, DT)
connected = full .- (QL .* QR)
ts = collect(0:M) .* DT
Emap = isempty(Emaps) ? zeros(0,0) : Matrix(reduce(hcat, Emaps)')
Jmap = isempty(Jmaps) ? zeros(0,0) : Matrix(reduce(hcat, Jmaps)')
jldsave(donef; theta_over_pi=theta_over_pi, t2p=t2p, N=N, ag=AG, mg=MG, center=CENTER,
        wilson=collect(wilson_endpoints()), R=R, bonds=[bL,bR], T=T, dt=DT, maxbond=mb,
        t=ts, full=full, connected=connected, QL=QL, QR=QR,
        heat_t=heat_t, Emap=Emap, Jmap=Jmap)
@printf("[%s] DONE  full(T)=%.4e%+.4eim  connected(T)=%.4e  -> %s\n",
        tag, real(full[end]), imag(full[end]), real(connected[end]), donef)

# =============================================================================
# Example (iii): charge-transfer detector correlator inside a Wilson-line state.
#
# Prepare |ψ⟩ = W|vac⟩ (a Wilson-line quench — NOT an energy eigenstate) and measure the
# correlator of the time-integrated current ("charge-transfer") detectors
#
#     Q(x) = ∫₀ᵀ j¹(x, t) dt ,     C(R, T) = ⟨ψ| Q(c−R) Q(c+R) |ψ⟩ ,
#
# for detector half-separations R. This is the genuinely two-time object
# ∫₀ᵀ∫₀ᵀ dt₁dt₂ ⟨ψ| j¹(c−R,t₁) j¹(c+R,t₂) |ψ⟩; the two evolutions are real many-body
# evolutions (no scalar-phase shortcut applies, since |ψ⟩ is not an eigenstate).
#
# Demonstrates `charge_transfer_correlator`:
#   1. a method cross-check (`:both`) at low resolution,
#   2. the cumulative curve C(R, t_m) vs upper limit t_m (free convergence scan in T),
#   3. a maxbond convergence check.
#
# Run from the project root:
#   SMOKE=1 julia --project=examples examples/wilson_charge_transfer.jl
#   AG=0.2 MG=0.0 N=64 NL=8 T=10 RLIST=6,12 NSTEPS=80 julia --project=examples examples/wilson_charge_transfer.jl
# =============================================================================
using Schwinger, SchwingerDetectors, JLD2, Printf

const SMOKE = get(ENV, "SMOKE", "0") == "1"
AG   = parse(Float64, get(ENV, "AG", "0.2"))
MG   = parse(Float64, get(ENV, "MG", "0.0"))
N    = parse(Int,     get(ENV, "N",  SMOKE ? "16" : "64"))
NL   = parse(Int,     get(ENV, "NL", SMOKE ? "4"  : "8"))       # Wilson-line length (sites)
T    = parse(Float64, get(ENV, "T",  SMOKE ? "3.0" : "10.0"))
NST  = parse(Int,     get(ENV, "NSTEPS", SMOKE ? "12" : "80"))  # grid intervals M
SUB  = parse(Int,     get(ENV, "SUBSTEPS", SMOKE ? "2" : "2"))  # TDVP substeps per gap
MAXB = parse(Int,     get(ENV, "MAXBOND", SMOKE ? "48" : "200"))
Rlist = [parse(Int, x) for x in split(get(ENV, "RLIST", SMOKE ? "3" : "6,12"), ",")]
OUT  = get(ENV, "OUT", joinpath(@__DIR__, "..", "data", "toolkit_examples"))
mkpath(OUT)

# 1. Vacuum → Wilson-line quench |ψ⟩ = W|vac⟩.
m  = build_model(; N = N, ag = AG, mg = MG)
c  = N ÷ 2
l1 = c - NL ÷ 2; l2 = l1 + NL
wl = wilson_line_quench(m, l1, l2; gs_kwargs = (; energy_tol = 1e-9))
psi = wl.state
@printf("wilson line on sites %d…%d;  charge-transfer detectors Q(c±R), R = %s\n", l1, l2, Rlist)

# 2. Method cross-check at low resolution (O(M²) matrix vs O(M) accumulated).
R0 = first(Rlist)
chk = charge_transfer_correlator(psi, min(T, 3.0), (c - R0, c + R0);
                                 nsteps = min(NST, 12), substeps = SUB,
                                 two_site = true, maxbond = MAXB, method = :both)
@printf("cross-check (R=%d): max|accumulated − matrix| = %.3e\n", R0, chk.discrepancy)

# 3. Production: cumulative curve C(R, t_m) for each R (accumulated, O(M)).
tgrid = Float64[]
full = Dict{Int,Vector{ComplexF64}}(); conn = Dict{Int,Vector{ComplexF64}}()
for R in Rlist
    res = charge_transfer_correlator(psi, T, (c - R, c + R);
                                     nsteps = NST, substeps = SUB, two_site = true,
                                     maxbond = MAXB, method = :accumulated)
    global tgrid = res.t
    full[R] = res.full; conn[R] = res.connected
    @printf("R=%2d:  C_full(T)=(% .4e,% .4e)   C_conn(T)=% .4e\n",
            R, real(res.full[end]), imag(res.full[end]), real(res.connected[end]))
end

# 4. Bond-dimension convergence at the largest R.
res_lo = charge_transfer_correlator(psi, T, (c - last(Rlist), c + last(Rlist));
                                    nsteps = NST, substeps = SUB, two_site = true,
                                    maxbond = MAXB ÷ 2, method = :accumulated)
@printf("maxbond convergence (R=%d): |C(T; D=%d) − C(T; D=%d)| = %.3e\n",
        last(Rlist), MAXB, MAXB ÷ 2, abs(full[last(Rlist)][end] - res_lo.full[end]))

jldsave(joinpath(OUT, "wilson_charge_transfer$(SMOKE ? "_smoke" : "").jld2");
        N = N, ag = AG, mg = MG, l1 = l1, l2 = l2, Rlist = Rlist, nsteps = NST,
        substeps = SUB, maxbond = MAXB, t = tgrid,
        full = Dict(string(R) => full[R] for R in Rlist),
        connected = Dict(string(R) => conn[R] for R in Rlist),
        cross_check_discrepancy = chk.discrepancy)
println("wrote $(joinpath(OUT, "wilson_charge_transfer$(SMOKE ? "_smoke" : "").jld2"))")

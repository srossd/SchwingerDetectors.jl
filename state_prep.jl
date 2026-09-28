# =============================================================================
# state_prep.jl — state preparation for real-time Schwinger-model simulations.
#
# Every routine here is built from EXPORTED Schwinger.jl calls: Lattice, Hamiltonian,
# groundstate, loweststates (+ the `vacuumof`/`quasiparticle` accessors), wavepacket,
# WilsonLine, FermionField, act, and the sudden-quench primitive `rehost`/`quench`.
#
# As of the 2026-09-28 Schwinger.jl release these are all first-class exported features
# (they used to require reaching into internals — see UPSTREAM_WISHLIST.md, now resolved).
# =============================================================================

# -----------------------------------------------------------------------------
# Model construction
# -----------------------------------------------------------------------------

"""
    build_model(; N, F=1, q=1, ag=0.2, mg=0.0, theta2pi=0.0,
                  backend=:MPSKit, flavor_sym=false) -> NamedTuple

Bundle a lattice + Hamiltonian and the scalar parameters that made them. Units are
`g = 1`, so `ag = a·g` is the spacing and `mg = m/g` the mass (matching CLAUDE.md).

`N` may be an integer (finite lattice) or `Inf` (infinite, for wavepackets/correlators).
`theta2pi` may be a scalar or a per-site vector (a θ-step profile — see `string_quench`).

Returns `(; lat, H, N, F, q, ag, mg, theta2pi, backend, flavor_sym)`; pass the whole
thing to the state-prep routines below.
"""
function build_model(; N, F::Int = 1, q::Int = 1, ag::Real = 0.2, mg = 0.0,
                       theta2pi = 0.0, backend = :MPSKit, flavor_sym::Bool = false)
    lat = Lattice(N; F = F, q = q, a = ag, m = mg, θ2π = theta2pi, flavor_sym = flavor_sym)
    H   = Hamiltonian(lat; backend = backend)
    return (; lat, H, N, F, q, ag, mg, theta2pi, backend, flavor_sym)
end

"Vacuum / ground state of the model. Forwards `energy_tol`, `cutoff`, `bonddim`, … to `groundstate`."
prepare_groundstate(m; kwargs...) = groundstate(m.H; kwargs...)

# The sudden-quench primitive (re-host a wavefunction under a new Hamiltonian) is now the
# exported `rehost` / `quench` in Schwinger.jl — used directly below.

# -----------------------------------------------------------------------------
# Quench 1: static probe-charge "string" quench (θ-step)
# -----------------------------------------------------------------------------

"""
    string_quench(m; length_links, center=nothing, dtheta=1.0,
                  gs=nothing, gs_kwargs=(;)) -> NamedTuple

Suddenly insert a pair of static ±q probe charges a distance `length_links` apart,
realised as a `+dtheta` step in the per-link `θ/2π` profile between them (the exact
Gauss-law field of a ±q defect pair, but on a plain N-site lattice so the per-bond
current/energy-current detectors still work — the genuine defect-*site* build inserts
inert sites that those operators cannot address).

Prepares the vacuum of `m` (or reuses a passed `gs`), builds the stepped Hamiltonian,
and returns the quench state (vacuum re-hosted under it). Finite lattice only.

Returns `(; state, lat, H, gs, links)` where `links = (l1, l2)` is the stepped range.
"""
function string_quench(m; length_links::Int, center = nothing, dtheta::Real = 1.0,
                        gs = nothing, gs_kwargs = (;))
    isfinite(m.lat.N) || throw(ArgumentError("string_quench needs a finite lattice"))
    N  = Int(m.lat.N)
    gs = isnothing(gs) ? groundstate(m.H; gs_kwargs...) : gs
    c  = isnothing(center) ? N ÷ 2 : center
    l1 = c - length_links ÷ 2
    l2 = l1 + length_links
    (2 ≤ l1 < l2 ≤ N) || throw(ArgumentError("string [$l1,$l2) out of 1..$N"))
    θ = m.theta2pi isa Number ? fill(Float64(m.theta2pi), N) : collect(Float64.(m.theta2pi))
    for n in l1:(l2 - 1); θ[n] += dtheta; end
    latq = Lattice(N; F = m.F, q = m.q, a = m.ag, m = m.mg, θ2π = θ, flavor_sym = m.flavor_sym)
    Hq   = Hamiltonian(latq; backend = m.backend)
    return (; state = rehost(gs, Hq), lat = latq, H = Hq, gs = gs, links = (l1, l2))
end

# -----------------------------------------------------------------------------
# Quench 2: Wilson-line quench (act a gauge-invariant χ†…χ string on the vacuum)
# -----------------------------------------------------------------------------

"""
    wilson_line_quench(m, start, finish; gs=nothing, flavor=1, conjugate=false,
                       normalize=true, gs_kwargs=(;)) -> NamedTuple

Create a charge pair joined by an electric-flux string by applying the spatial Wilson
line `WilsonLine(lat, conjugate, flavor, start, finish)` to the vacuum, then evolve
under the *unmodified* Hamiltonian. On a θ=π vacuum (`theta2pi=0.5`) this instead flips
the electric field over `[start,finish)` and creates a soliton–antisoliton pair;
`conjugate=true` selects the tensionless partner.

Returns `(; state, gs, op)`.
"""
function wilson_line_quench(m, start::Int, finish::Int; gs = nothing, flavor::Int = 1,
                            conjugate::Bool = false, normalize::Bool = true, gs_kwargs = (;))
    gs = isnothing(gs) ? groundstate(m.H; gs_kwargs...) : gs
    W  = WilsonLine(m.lat, conjugate, flavor, start, finish; backend = m.backend)
    st = act(W, gs)
    normalize && normalize!(st.psi)
    return (; state = st, gs = gs, op = W)
end

# -----------------------------------------------------------------------------
# Quench 3: generic local-operator quench
# -----------------------------------------------------------------------------

"""
    local_operator_quench(gs, op; normalize=false) -> MPSKitState

Apply any (exported) local operator `op` to the state `gs` and return `op|gs⟩`.
Handy `op` builders (all exported): `FermionField(lat, site; backend, dagger, flavor)`
(changes fermion number), `WilsonLine(...)`, `ChargeCurrent(lat, bond; backend)`,
`EnergyCurrent(lat, site; backend)`. Leave `normalize=false` to preserve the operator's
magnitude (needed for correlators); set `true` for a genuine normalised quench state.

NOTE: this works on a plain finite-lattice state. Applying a per-bond/site current
operator to a WindowMPS (infinite background) is NOT supported — see README.
"""
function local_operator_quench(gs::MPSKitState, op; normalize::Bool = false)
    _assert_not_window(gs, "local_operator_quench")
    st = act(op, gs)
    normalize && normalize!(st.psi)
    return st
end

# -----------------------------------------------------------------------------
# Quench 4: moving soliton(s) and mesons — momentum QP → localized wavepacket
# -----------------------------------------------------------------------------

"""
    aligned_support(center, len) -> UnitRange

A length-`len` site window centred on `center`, snapped to start on an odd site and end
on an even site (the unit-cell alignment `wavepacket` requires).
"""
function aligned_support(center::Int, len::Int)
    lo = center - len ÷ 2 + 1
    hi = center + len ÷ 2
    lo += iseven(lo) ? 1 : 0
    hi -= isodd(hi)  ? 1 : 0
    return lo:hi
end

"""
    moving_soliton(m; momentum, W, support, sigma=nothing, center=nothing,
                   window=true, dp=0.05, gs_kwargs=(;)) -> NamedTuple

A single soliton (the θ=π domain wall) launched at physical `momentum`, built as a
Gaussian wavepacket on a `WindowMPS` of `W` sites. Requires `m` to be an INFINITE
lattice at `theta2pi = 0.5`. Also returns the group velocity `vg` (from a ±`dp`
finite difference) so you can size the evolution window.

Returns `(; state, qp, E1, vg, vacuum)`.
"""
function moving_soliton(m; momentum::Real, W::Int, support,
                        sigma = nothing, center = nothing, window::Bool = true,
                        dp::Real = 0.05, gs_kwargs = (;))
    isinf(m.lat.N) || throw(ArgumentError("moving_soliton needs an infinite lattice (N=Inf)"))
    res = loweststates(m.H, 2; solitons = true,
                       momentum = [momentum, momentum + dp, momentum - dp], gs_kwargs...)
    v1      = vacuumof(res)                                                  # accessor: hides nesting
    soliton = quasiparticle(res, 2; momentum = 1, kind = :soliton)          # soliton @ +p
    E1 = real(energy(soliton))
    vg = (real(energy(quasiparticle(res, 2; momentum = 2, kind = :soliton))) -
          real(energy(quasiparticle(res, 2; momentum = 3, kind = :soliton)))) / (2dp)
    wp = wavepacket([soliton], W; supports = [support], sigmas = [sigma],
                    centers = [center], gauge = :symmetric, window = window)
    return (; state = wp, qp = soliton, E1 = E1, vg = vg, vacuum = v1)
end

"""
    moving_meson(m; momentum, W, support, sigma=nothing, center=nothing,
                 window=true, dp=0.05, gs_kwargs=(;)) -> NamedTuple

A single neutral quasiparticle (the Schwinger boson / meson) launched at `momentum`.
Same machinery as [`moving_soliton`] but WITHOUT `solitons=true`; requires an INFINITE
lattice (typically `theta2pi = 0.0`). Returns `(; state, qp, E1, vg, vacuum)`.
"""
function moving_meson(m; momentum::Real, W::Int, support,
                      sigma = nothing, center = nothing, window::Bool = true,
                      dp::Real = 0.05, gs_kwargs = (;))
    isinf(m.lat.N) || throw(ArgumentError("moving_meson needs an infinite lattice (N=Inf)"))
    res = loweststates(m.H, 2; momentum = [momentum, momentum + dp, momentum - dp], gs_kwargs...)
    vac   = vacuumof(res)
    meson = quasiparticle(res, 2; momentum = 1)              # accessor: hides momentum-list nesting
    E1 = real(energy(meson))
    vg = (real(energy(quasiparticle(res, 2; momentum = 2))) -
          real(energy(quasiparticle(res, 2; momentum = 3)))) / (2dp)
    wp = wavepacket([meson], W; supports = [support], sigmas = [sigma],
                    centers = [center], gauge = :symmetric, window = window)
    return (; state = wp, qp = meson, E1 = E1, vg = vg, vacuum = vac)
end

"""
    moving_wavepackets(qps, W; supports, sigmas=nothing, centers=nothing,
                       window=true) -> MPSKitState

Several quasiparticles (solitons and/or mesons) as one multi-packet `WindowMPS` — e.g.
a soliton–antisoliton pair, or two mesons set up to collide. `qps` is a vector of the
`qp` fields returned by [`moving_soliton`]/[`moving_meson`]; `supports` a vector of
disjoint site ranges (use [`aligned_support`]). Thin wrapper over the exported
multi-state `wavepacket`.
"""
function moving_wavepackets(qps::AbstractVector, W::Int; supports,
                            sigmas = nothing, centers = nothing, window::Bool = true)
    n = length(qps)
    sig = isnothing(sigmas)  ? fill(nothing, n) : sigmas
    cen = isnothing(centers) ? fill(nothing, n) : centers
    return wavepacket(collect(qps), W; supports = collect(supports),
                      sigmas = collect(sig), centers = collect(cen),
                      gauge = :symmetric, window = window)
end

# -----------------------------------------------------------------------------
# internal helper
# -----------------------------------------------------------------------------

# A state hosted on an infinite lattice (a WindowMPS wavepacket, or a finite window on an
# infinite background) cannot take a per-bond/site operator via act/expectation — the
# operator is built on the infinite lattice. Only plain finite-lattice states can.
_iswindow(st::MPSKitState) = isinf(lattice(st).N)
function _assert_not_window(st::MPSKitState, who)
    _iswindow(st) && throw(ArgumentError(
        "$who: cannot apply a per-bond/site operator to an infinite-lattice (WindowMPS / " *
        "wavepacket) state. Read densities with the vectorized functions instead " *
        "(chargecurrents/energycurrents/energy_densities/…) — see README."))
    return nothing
end

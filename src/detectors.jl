# =============================================================================
# detectors.jl — detector-operator 1-point and 2-point functions.
#
# The "detectors" are the conserved currents (vector charge current j¹ = ChargeCurrent,
# energy current 𝒥 = T⁰¹ = EnergyCurrent) and the densities (energy, ρ, E, pseudoscalar,
# momentum). A detector sits at a FIXED lattice location and reads its operator vs time.
#
#   • 1-point ⟨O(x_D, t)⟩ — via the vectorized densities (works on WindowMPS AND finite
#     lattices), or an operator/`expectation` variant for finite lattices.
#   • 2-point ⟨O(x_D, t) O'(x_S, 0)⟩ — now the built-in `correlator2pt`; `detector_2pt`
#     here is a thin, named wrapper over it.
# =============================================================================

# -----------------------------------------------------------------------------
# 1-point functions
# -----------------------------------------------------------------------------

"""
    detector_1pt(state, T, locations; kind=:energycurrent, nsteps, kwargs...)
        -> (times, Dict(location => series))

Evolve `state` to `T` and return the time series of density `kind` at each fixed
`locations` (site indices, or bond indices for `:current`/`:momentum`). Universal: works
on finite-lattice states AND WindowMPS wavepackets (it reads the vectorized density —
`kind` is a built-in `standard_densities` menu name). Extra `kwargs` forward to
[`evolve_and_measure`] (`two_site`, `maxbond`, `subtract`, `save`, …).
"""
function detector_1pt(state, T::Real, locations; kind::Symbol = :energycurrent,
                      nsteps::Int, kwargs...)
    _, obs = evolve_and_measure(state, T; nsteps = nsteps, densities = [kind], kwargs...)
    M = density_map(obs, String(kind))                 # nt × nx
    series = Dict(x => M[:, x] for x in locations)
    return times(obs), series
end

"""
    detector_1pt_operator(state, T, sites; kind=:energycurrent, nsteps,
                          two_site=true, maxbond=300) -> (times, Matrix)

Finite-lattice detector via the per-site/bond conserved-current OPERATOR and
`expectation` (the `step4c`-style measurement): builds `EnergyCurrent`/`ChargeCurrent`/
`MomentumDensity` at each detector site and evaluates it every step. Returns
`(times, values)` with `values[t_index, det_index]`.

Only valid for a plain finite-lattice state (a per-bond operator cannot be applied to a
WindowMPS). For whole-lattice profiles use the vectorized `energycurrents(ψ; pad=true)` /
`detector_1pt` instead. Errors clearly on a WindowMPS.
"""
function detector_1pt_operator(state::MPSKitState, T::Real, sites; kind::Symbol = :energycurrent,
                               nsteps::Int, two_site::Bool = true, maxbond = 300)
    _assert_not_window(state, "detector_1pt_operator")
    lat = state.hamiltonian.lattice
    opof = if kind === :energycurrent
        s -> EnergyCurrent(lat, s; backend = :MPSKit)
    elseif kind === :current
        s -> ChargeCurrent(lat, s; backend = :MPSKit)
    elseif kind === :momentum
        s -> MomentumDensity(lat, s; backend = :MPSKit)
    else
        throw(ArgumentError("detector_1pt_operator kind must be :energycurrent, :current, or :momentum"))
    end
    ops  = Dict(s => opof(s) for s in sites)
    obsf = (ψ, t) -> [real(expectation(ops[s], ψ)) for s in sites]
    mb   = two_site ? maxbond : nothing
    _, obs = evolve(state, T; nsteps = nsteps, two_site = two_site, maxlinkdim = mb,
                    observable = Dict("det" => obsf))
    return times(obs), permutedims(reduce(hcat, collect(obs.det)))   # nt × ndet
end

# -----------------------------------------------------------------------------
# 2-point functions
# -----------------------------------------------------------------------------

"""
    detector_2pt(vac, source_op, detector_op, tgrid; connected=true, kwargs...)
        -> (tgrid, C)

Real-time two-time correlator `C(t) = ⟨vac| detector_op(t) source_op(0) |vac⟩` on a
FINITE lattice, where `vac` is an energy eigenstate (the vacuum). A thin named wrapper
over the built-in `correlator2pt(vac, detector_op, source_op, tgrid; …)`, which assembles
it from Schwinger.jl's own `act`/`evolve`/`dot` and removes the trivial vacuum phase
`e^{iE₀t}` with `E₀ = ⟨vac|H|vac⟩`. `detector_op` is assumed Hermitian (pass its dagger
otherwise). Extra `kwargs` (`two_site`, `maxlinkdim`, `E0`, `nsteps`) forward to
`correlator2pt`. `tgrid` must be sorted, non-negative.

Operators are any exported locals, e.g. `ChargeCurrent(lat, bond; backend=:MPSKit)`,
`EnergyCurrent(lat, site; backend=:MPSKit)`, `WilsonLine(...)`.
"""
function detector_2pt(vac::SchwingerState, source_op::SchwingerOperator,
                      detector_op::SchwingerOperator, tgrid; connected::Bool = true, kwargs...)
    C = correlator2pt(vac, detector_op, source_op, tgrid; connected = connected, kwargs...)
    return collect(tgrid), C
end

"""
    detector_2pt_equal_time(state, bonds; kind=:current, connected=true)
        -> Dict((bR,bL) => value)

Equal-time connected two-current correlator on the CURRENT state (e.g. a snapshot of a
quench), `⟨state| O(bR) O(bL) |state⟩ − ⟨O(bR)⟩⟨O(bL)⟩`, using the exported operator
product. Finite lattice only. `bonds` is a list of `(bR, bL)` pairs; `kind` ∈
(`:current`, `:energycurrent`). This is the correct object for a NON-eigenstate quench
state (where the single-phase `correlator2pt` trick does not apply).
"""
function detector_2pt_equal_time(state::MPSKitState, bonds; kind::Symbol = :current,
                                 connected::Bool = true)
    _assert_not_window(state, "detector_2pt_equal_time")
    lat = state.hamiltonian.lattice
    op  = kind === :current       ? (s -> ChargeCurrent(lat, s; backend = :MPSKit)) :
          kind === :energycurrent ? (s -> EnergyCurrent(lat, s; backend = :MPSKit)) :
          throw(ArgumentError("kind must be :current or :energycurrent"))
    one = kind === :current ? chargecurrents(state) : energycurrents(state)
    out = Dict{Tuple{Int,Int},Float64}()
    for (bR, bL) in bonds
        full = real(expectation(op(bR) * op(bL), state))
        out[(bR, bL)] = connected ? full - real(one[bR]) * real(one[bL]) : full
    end
    return out
end

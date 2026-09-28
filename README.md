# Schwinger toolkit — consolidated real-time capabilities

A small, documented library that collects the recurring capabilities scattered across
`Schwinger/`, `Schwinger/correlators/`, and `confinement/` into **three short files**,
built from **exported `Schwinger.jl` calls**. Load it with

```julia
include("Schwinger/toolkit/SchwingerToolkit.jl")
using .SchwingerToolkit
```

and run the wiring check with `SMOKE=1 julia --project=. Schwinger/toolkit/demo.jl`.

Full API docs (Documenter.jl, mirroring `Schwinger.jl`) live in `docs/`; build with
`julia --project=docs docs/make.jl`.

Units are `g = 1` throughout: `ag = a·g`, `mg = m/g` (matching `CLAUDE.md`).

> **As of the 2026-09-28 Schwinger.jl release**, the package absorbed the whole wishlist
> (`UPSTREAM_WISHLIST.md`): the sudden-quench primitive (`rehost`/`quench`), the
> `standard_densities` observable menu, the real-time two-point correlator
> (`correlator2pt`), local operator insertion on windows (`apply_local`), the
> `vacuumof`/`quasiparticle` accessors, `dispersion`/`groupvelocity`, `chargeprofile`,
> `normalize`, `savestate`/`loadstate`, and `energy_densities(…; convention=:bond)`. The
> toolkit now just calls these; it is thin convenience glue, not a workaround layer.

## The three files

| file | capability | key functions |
|------|------------|---------------|
| `state_prep.jl` | prepare initial states | `build_model`, `prepare_groundstate`, `string_quench`, `wilson_line_quench`, `local_operator_quench`, `moving_soliton`, `moving_meson`, `moving_wavepackets` |
| `measure.jl` | evolve + record densities | `evolve_and_measure`, `times`, `density_map`, `vacuum_tile` |
| `detectors.jl` | detector 1-pt / 2-pt | `detector_1pt`, `detector_1pt_operator`, `detector_2pt`, `detector_2pt_equal_time` |

`SchwingerToolkit.jl` is the module that includes and re-exports them; `demo.jl` is a
tiny end-to-end example. Two focused examples live in `examples/`:

- `examples/energy_1pt_probe_charges.jl` — **(i)** energy one-point function ⟨𝒥(x_D,t)⟩
  after inserting a static probe-charge pair.
- `examples/charge_charge_wilson.jl` — **(ii)** connected equal-time charge–charge
  correlator ⟨ψ(t)| j¹(c+R) j¹(c−R) |ψ(t)⟩_c *inside* the Wilson-line-quenched state
  |ψ⟩ = W|vac⟩, measured on each snapshot with `detector_2pt_equal_time`.

## 1. State preparation (`state_prep.jl`)

Start with a model, then pick a preparation. Everything here is clean exported API.

```julia
m = build_model(; N = 128, ag = 0.2, mg = 0.0)          # finite; N = Inf for wavepackets
gs = prepare_groundstate(m; energy_tol = 1e-9)          # = groundstate(m.H; …)
```

- **Ground state** — `prepare_groundstate(m; …)`.
- **Static probe-charge "string" quench** — `string_quench(m; length_links, dtheta=1)`.
  Inserts a ±q pair as a `θ/2π` step (Gauss-law field of a defect pair) on a *plain*
  N-site lattice, so the per-bond current/energy-current detectors still work. Returns
  `(; state, lat, H, gs, links)`.
- **Wilson-line quench** — `wilson_line_quench(m, start, finish; conjugate, flavor)`.
  Applies the gauge-invariant `χ†…χ` string `WilsonLine` to the vacuum (a meson at θ=0;
  a soliton–antisoliton bubble at θ=π, where `conjugate=true` picks the tensionless one).
- **Generic local-operator quench** — `local_operator_quench(gs, op)` for any exported
  local `op` (`FermionField`, `WilsonLine`, `ChargeCurrent`, …).
- **Moving soliton** — `moving_soliton(m; momentum, W, support, …)` (needs `N=Inf`,
  `theta2pi=0.5`). Builds the momentum-eigenstate soliton and Gaussian-envelopes it into
  a `WindowMPS` wavepacket; returns the group velocity `vg` too.
- **Moving meson** — `moving_meson(m; momentum, W, support, …)` (needs `N=Inf`), the
  neutral-quasiparticle analogue.
- **Multiple / colliding packets** — `moving_wavepackets([qpA, qpB], W; supports=[…])`.

Helper: `aligned_support(center, len)` builds a unit-cell-aligned site window.

## 2. Evolve while measuring densities (`measure.jl`)

```julia
evolved, obs = evolve_and_measure(state, T; nsteps = 100, two_site = true, maxbond = 128,
                                  densities = [:energy, :electricfield, :current,
                                               :energycurrent, :pseudoscalar])
t   = times(obs)                       # Vector of times
Jxt = density_map(obs, :energycurrent) # nt × nx spacetime array
```

The density menu is the built-in `standard_densities` (all exported observables,
WindowMPS-safe): `:charge` (ρ = j⁰), `:occupation`, `:electricfield`, `:energy`,
`:current` (j¹ per bond), `:energycurrent` (𝒥 = T⁰¹ per site), `:pseudoscalar`,
`:momentum`. (For the 𝒥-consistent energy density call
`energy_densities(state; convention=:bond)` directly.)

- `subtract = Dict(:energy => vacuum_tile(vac, :energy, W), …)` stores only the excess
  over an infinite-lattice vacuum (use `vacuum_tile` to tile the staggered background).
- `save = (state, obs) -> jldsave(...)`, `save_every = 10` for checkpointing.
- `grow = true` enables adaptive window growth for a `WindowMPS`.

## 3. Detector 1-point / 2-point functions (`detectors.jl`)

**1-point** ⟨O(x_D, t)⟩ at fixed detector locations:

```julia
t, series = detector_1pt(state, T, [xL, xR]; kind = :energycurrent, nsteps = 100,
                         two_site = true, maxbond = 128)   # universal (works on WindowMPS)
```

`detector_1pt` reads the vectorized density (works everywhere). `detector_1pt_operator`
is the operator/`expectation` variant (finite lattices only) matching the historical
`step4c_detector.jl`.

**2-point** ⟨O(x_D, t) O'(x_S, 0)⟩ on a **finite lattice** — this is the built-in
`correlator2pt` (Schwinger.jl assembles it from its own `act`/`evolve`/`dot` and removes
the vacuum phase `e^{iE₀t}` with `E₀ = ⟨vac|H|vac⟩`). `detector_2pt` is a thin named
wrapper:

```julia
lat = m.lat
# either the built-in directly …
C = correlator2pt(vac, ChargeCurrent(lat, bD), ChargeCurrent(lat, bS),
                  collect(0.0:0.1:10.0); connected = true, two_site = true, maxlinkdim = 128)
# … or the toolkit wrapper (source, detector) order + returns (tgrid, C):
t, C = detector_2pt(vac,
                    ChargeCurrent(lat, bS; backend = :MPSKit),   # source at bond bS (t=0)
                    ChargeCurrent(lat, bD; backend = :MPSKit),   # detector at bond bD (t)
                    collect(0.0:0.1:10.0); connected = true, two_site = true, maxlinkdim = 128)
```

`correlator2pt` assumes an **eigenstate** reference (the vacuum). For the equal-time
⟨j¹ j¹⟩ *inside* a non-eigenstate quench state as it evolves, use
`detector_2pt_equal_time(state, [(bR,bL), …])` on each snapshot.

---

## Caveats and remaining edges

Most of the old workarounds are gone (the wishlist was implemented upstream). What's left:

1. **Per-bond/site current operators on an infinite-lattice (WindowMPS / wavepacket)
   state.** `expectation(EnergyCurrent(lat, n), ψ)` / `act(ChargeCurrent(lat, b), ψ)` only
   work when `ψ` is a *plain finite-lattice* state. On a `WindowMPS` (soliton / meson /
   window-vacuum) use the vectorized densities (`chargecurrents`, `energycurrents(…;
   pad=true)`, `energy_densities`, …) — `detector_1pt` and `evolve_and_measure` do this,
   so they are safe on both; the operator routines (`detector_1pt_operator`,
   `detector_2pt`, `detector_2pt_equal_time`) error clearly on a window state.
   *Now partly liftable:* single-site operators (charge density, phases) **can** be applied
   to a window with the new `apply_local`; per-*bond* bilinears (a current) still can't.

2. **Real-time correlators need an eigenstate reference.** `correlator2pt` / `detector_2pt`
   remove the vacuum phase `e^{iE₀t}` with the single number `E₀ = ⟨vac|H|vac⟩`, which is
   exact only when the reference is an energy eigenstate. For a two-time correlator inside a
   *quench* state, use `detector_2pt_equal_time` on snapshots. The genuine infinite-vacuum
   `⟨0|jᵘ(t,x)…|0⟩` n-point *tabulations* over a box still live in
   `Schwinger/correlators/` (they need per-bond insertions on a window — see item 1).

3. **Diagonal phase-unitary source** `W_f = ∏ₖ e^{iθₖ Qₖ}` (anomaly detector) is still a
   hand-built `FiniteMPO`; but with `apply_local` you can now apply the single-site phase
   factors to a window directly.

4. **`DefectCharge`-*site* quench** (`insert_defect`) is exported and clean, but its inert
   sites break per-bond detectors — hence `string_quench` uses the equivalent θ-step.

5. **Adjoint QCD₂** (`adjoint/`) is a *different package* (`AdjointQCD2.jl`); this toolkit
   is Schwinger-only.

## Provenance

Consolidated from: `step4c_detector.jl`, `soliton_charge_detection.jl`,
`meson_axial_ward.jl`, `ward_sim.jl`, `ward_wilson.jl`, `confinement/string_quench.jl`,
`correlators/{correlators,run_2pt,soliton_wilson,soliton_corr}.jl`, and the
`Schwinger.jl` unified API (`states.jl`, `currents.jl`, `timeevolution.jl`, `api.jl`).

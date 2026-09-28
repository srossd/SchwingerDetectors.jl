# Schwinger.jl wishlist — missing features, papercuts, and suspected bugs

> ## ✅ RESOLVED — Schwinger.jl release 2026-09-28
> The developer implemented essentially this entire list (commits `9d21602`, `c72cb76`,
> `1de6585`). New exported API:
> - **§2.1** `apply_local(state, op, site)` — single-site operator insertion on a window
>   (multi-site bond bilinears still pending).
> - **§2.2** `correlator2pt` / `correlator2pt_states` — real-time two-point functions.
> - **§2.3** generic `chargecurrents`/`energycurrents`/`momentumdensities`/`totalmomentum`
>   for ED & ITensors.  **§2.5** momentum density unified + exported.
> - **§2.7** `dispersion(H, momenta)` / `groupvelocity(H, p)`.  **§2.8** `normalize`(`!`),
>   `savestate`/`loadstate`.  **§2.9** `rehost`/`quench`.  **§2.10** `standard_densities`.
> - **§1.1** `energy_densities(…; convention=:bond)` (the 𝒥 partner).  **§1.2/1.3**
>   pseudoscalar docstrings + a regression test locking cross-backend agreement (no bug).
> - **§3.1** `vacuumof`/`quasiparticle` accessors.  **§3.3** `chargeprofile`.
>   **§3.4** `pad=true` on `energycurrents`/`momentumdensities`.
>
> The toolkit in this folder has been updated to call these built-ins. The original
> analysis is kept below for provenance.

Compiled while consolidating the real-time detector-operator work into
`Schwinger/toolkit/`. This is a punch list for **upstream `Schwinger.jl`**
(`~/.julia/dev/Schwinger`): things that, if they existed, would have let the toolkit (and
the `Schwinger/`, `Schwinger/correlators/`, `confinement/` scripts) avoid reaching past
the exported API.

**Confidence legend** — each item is tagged:
- 🐞 **BUG / correctness footgun** — a wrong or surprising result is possible. I have *not*
  written a failing test for most of these; treat as "verify then fix". Confidence noted inline.
- 🧩 **MISSING FEATURE** — capability that only exists as hand-rolled internals downstream.
- ✂️ **PAPERCUT** — works, but the ergonomics force boilerplate or invite mistakes.

Line references are `file:line` in the dev checkout as of this writing.

---

## 1. Correctness / suspected bugs

### 1.1 🐞 `energy_densities` is not the density that satisfies the 𝒥 continuity equation
*Confidence: medium-high — documented upstream, not re-derived here.*

The built-in `energy_densities` uses a **site-centered** splitting, whereas the exact
lattice continuity `∂_t h_n = 𝒥_n − 𝒥_{n+1}` holds for a **bond-centered** density
`hₙ = ⟨GK(n)⟩ + ⟨Hop(n)⟩ + ½⟨M(n)⟩ + ½⟨M(n+1)⟩`. The two differ by a total derivative, so
pointwise energy-continuity checks fail with the built-in and the scripts rebuild the
consistent density by hand. This is stated in the upstream examples README
(`examples/detector_operators/README.md:76-77`).

**Ask:** either export the bond-consistent density used by `EnergyCurrent` (e.g.
`energy_densities(state; convention=:bond)`), or document loudly that `energy_densities`
is *not* the partner of `EnergyCurrent`. As-is it's a trap for anyone doing local energy
conservation.

### 1.2 🐞 `pseudoscalardensity` naming vs the naive `j¹ ∝ (−1)ⁿ P` shortcut
*Confidence: medium — documented, but a naming/semantics trap.*

`pseudoscalardensity` is the **site-symmetrized staggered** pseudoscalar, so the textbook
shortcut `j¹ = κ(−1)ⁿ · pseudoscalardensity` does **not** reproduce the current (it gives a
bond difference). Anyone porting continuum formulas will get a wrong current and a
confusing factor. **Ask:** document the exact operator each `*density` returns, and/or
provide the un-symmetrized bond bilinear as a separate exported observable (see 2.4).

### 1.3 🐞 `pseudoscalardensities` has two methods; confirm they agree
*Confidence: low — flagging for a test, not a demonstrated failure.*

There is a generic `pseudoscalardensities(state::SchwingerState)` (`states.jl:1551`) **and**
an `MPSKitState` override (`currents.jl:223`). If these use different conventions
(site-symmetrized vs bond) an ED↔MPSKit comparison would silently disagree. **Ask:** a
cross-backend regression test asserting they return the same thing on a shared state.

---

## 2. Missing features (only exist as downstream internals)

### 2.1 🧩 Local operators cannot be applied to a `WindowMPS`
*The single biggest gap — it's the reason `correlators/correlators.jl` exists.*

A global `FiniteMPO * WindowMPS` product is unsupported, so **no** per-bond/site operator
(`ChargeCurrent`, `EnergyCurrent`, `MomentumDensity`, `WilsonLine`, `FermionField`) can be
`act`-ed or `expectation`-ed on a wavepacket / window-vacuum state. Downstream this forces
two workarounds:
- read observables only through the vectorized `chargecurrents`/`energycurrents`/… (fine,
  but read-only — you can't build `O|ψ⟩`);
- for genuine correlators, hand-build the operator as a TensorKit `TensorMap` and gate it
  onto the window with `MPSKit.svd_trunc!` / `MPSKit._transpose_front`, manipulating
  `ψ.window`, `ψ.AC`, `ψ.AR`, `ψ.left_gs`, `ψ.right_gs` directly
  (`correlators/correlators.jl:87-180`).

**Ask:** an exported `act(op, ::MPSKitState)` / `apply_local(state, op, site)` that works on
a `WindowMPS` for the local operators, applying the operator on the window tensors and
leaving the infinite wings untouched. This alone would delete ~200 lines of internal
tensor surgery from the correlator suite.

### 2.2 🧩 No exported real-time n-point correlator machinery
Every correlator script re-implements the same pipeline: window vacuum → insert local `O`
→ TDVP-evolve → insert → overlap → multiply by a **self-calibrated** vacuum phase
`e^{iE₀t}` obtained from `⟨0|e^{−iHt}|0⟩` (`correlators/correlators.jl:74-82`,
`run_2pt.jl`, `run_3pt_task.jl`). The E₀ calibration exists only because the extensive
window energy isn't otherwise available accurately.

**Ask:** export (a) the vacuum-phase / E₀ calibration as a helper, and (b) a
`correlator(vac, insertions; tgrid, …)` driver for 2-/3-point functions. The toolkit's
`detector_2pt` covers the *finite-lattice* two-time case with `act`/`evolve`/`dot`, but the
*infinite-background* tabulations can't be done without the internals above.

### 2.3 🧩 Vectorized current observables are MPSKit-only (cross-backend gap)
*Confidence: high — verified by grep.*

`chargecurrents` (`currents.jl:259`), `energycurrents` (`currents.jl:286`), and
`momentumdensities` (`currents.jl:424`) have **only** an `::MPSKitState` method. ED and
ITensors states have no vectorized whole-lattice current/energy-current/momentum profile
(you must loop `expectation(ChargeCurrent(lat,b), state)` by hand). `pseudoscalardensities`
*does* have a generic fallback, so the coverage is inconsistent even within the family.

**Ask:** generic `::SchwingerState` fallbacks (loop the per-bond operator) so
`energycurrents`/`chargecurrents`/`momentumdensities` work on every backend — important for
ED cross-checks.

### 2.4 🧩 No exported un-symmetrized bond bilinears (density/current/pseudoscalar on a bond)
The correlator engine needs the **charge density on a site**, the **current on a bond**,
and the **pseudoscalar bilinear on a bond** as standalone local operators; it hand-writes
all three (`correlators/correlators.jl:87-118`, comps 0/1/2). **Ask:** export these local
one-/two-site operators (they're the natural companions of `ChargeCurrent`).

### 2.5 🧩 No exported momentum-distribution moments
`totalmomentum(wp)` gives the mean, but ⟨P²⟩ / Var(P) require a hand-built total-P MPO from
`Schwinger.get_mpskit_spaces`, `Schwinger._momdens_prefactor`, `FiniteMPOHamiltonian`
(`examples/detector_operators/momentum_density_moments.jl`). **Ask:** export
`momentmomentum(wp, k)` or a `momentum_distribution(wp)`.

### 2.6 🧩 No exported diagonal phase-unitary source `∏ₖ e^{iθₖ Qₖ}`
The anomaly detector builds it as a raw `Schwinger.MPSKitOperator(lat, FiniteMPO(...),
universe)` (`correlators/detector_anomaly.jl:81-96`). **Ask:** export a
`ChargePhase(lat, θ)` / `gauge_rotation` operator; it's a natural quench source.

### 2.7 🧩 No exported group-velocity / dispersion helper
Every wavepacket script computes `v_g = (E(p+dp) − E(p−dp)) / 2dp` by calling
`loweststates` at three momenta and finite-differencing by hand (soliton/meson scripts,
`step4a_sumrule_dispersion.jl`). **Ask:** `groupvelocity(H, p; dp=…)` and/or a
`dispersion(H, ps)` band helper.

### 2.8 🧩 No exported state-level `normalize` / `save` / `load`
- `act` returns an **unnormalized** state and there is no `normalize(state)` /
  `normalize!(state)` — downstream reaches into `normalize!(state.psi)` everywhere.
- Persisting a state means `jldsave(...; psi = state.psi)` then rebuilding by hand as
  `MPSKitState(H, load(f,"psi"), defects)` on reload (`confinement/string_quench.jl:83`,
  `correlators/soliton_wilson.jl:45`). There's no `savestate`/`loadstate` that round-trips
  the Hamiltonian/defect metadata.

**Ask:** `LinearAlgebra.normalize(!)` methods for `SchwingerState`, and
`savestate(path, state)` / `loadstate(path, H)`.

### 2.9 🧩 No exported "quench" primitive (re-host a state under a new Hamiltonian)
The core sudden-quench move — take `|ψ⟩` from one Hamiltonian and evolve it under another —
is written as the low-level `MPSKitState(Hnew, ψ.psi, ψ.defects)` everywhere, and
`step2_axial_ward.jl:19-21` even hand-writes a `_withH` that dispatches over
`EDState`/`ITensorState`/`MPSKitState` (reaching into `.coeffs`/`.psi`/`.net_charge`).
**Ask:** an exported `rehost(state, Hnew)` / `quench(state, Hnew)` with methods for all
three backends (the toolkit ships a one-backend `rehost`; upstream should own the
cross-backend one).

### 2.10 🧩 No built-in observable menu for `evolve`
`evolve(...; observable=Dict(...))` is great, but everyone hand-writes the same
`"charge" => (ψ,t)->real.(vec(charges(ψ)))`, `"current" => (ψ,t)->chargecurrents(ψ)`, …
dictionary. **Ask:** a helper like `standard_densities([:charge,:current,:energy,…])` that
returns the callback dict (this is exactly what the toolkit's `DENSITY_FUNCS` /
`evolve_and_measure` provide — a good candidate to upstream).

---

## 3. API ergonomics / papercuts

### 3.1 ✂️ `loweststates` return shape differs for solitons vs ordinary QPs
For `solitons=true`, a quasiparticle is `res[2][j][1]` (extra nesting); for an ordinary
band it's `res[2][j]`. Indexing the wrong depth is a silent footgun (you get a container,
not a state). **Ask:** a uniform return type, or accessor helpers
(`quasiparticle(res, j)`, `vacuumof(res)`).

### 3.2 ✂️ Three coexisting backend-selection styles
`Hamiltonian(lat; backend=:MPSKit)` (symbol), `Hamiltonian(lat, MPSKitBackend())`
(positional object), and `set_default_backend(:MPSKit)` (global) all appear across the
scripts, plus ad-hoc `backend isa EDBackend ? … : …` branching for backend-agnostic code.
**Ask:** pick one blessed style in the docs; consider a `with_backend(:MPSKit) do … end`
scope for agnostic drivers.

### 3.3 ✂️ `charges` (and friends) return a matrix needing `vec`
Multi-flavor shape leaks into F=1 code as `real.(vec(charges(ψ)))` boilerplate. **Ask:** a
1-D return for F=1, or a documented `chargeprofile(ψ)` that always returns a vector.

### 3.4 ✂️ Boundary-invalid detector operators throw instead of returning NaN
`EnergyCurrent` is valid for `2 ≤ n ≤ N−1`, `MomentumDensity` for `1 ≤ n ≤ N−2`; asking for
an edge site throws. For sweeping detectors across the whole lattice this forces manual
range guards. **Ask:** optionally return `NaN`/`missing` at invalid sites.

### 3.5 ✂️ `evolve` window growth is incompatible with defects
`evolve(...; grow=…)` throws when the state carries defects (`timeevolution.jl`). Reasonable
limitation, but it means the defect-site quench and adaptive windows are mutually exclusive
— worth documenting alongside the recommendation to use the θ-step string quench instead.

### 3.6 ✂️ Defect-*site* build silently disables per-bond detectors
`insert_defect` inserts inert sites the per-bond current/energy-current operators cannot
address, so any detector measurement on a defect-site state is wrong/impossible; the θ-step
(θ2π profile) equivalent is required instead. This is "working as intended" but is a sharp
edge with no guard. **Ask:** either make the per-bond operators defect-aware, or error when
a current operator is evaluated on a defect-site lattice.

### 3.7 ✂️ Cross-package `WilsonLine` signature drift (Schwinger vs AdjointQCD2)
Schwinger's `WilsonLine(lat, conj::Bool, flavor::Int, start, finish; backend)` vs
AdjointQCD2's `WilsonLine(lat, start, finish; fermionparity)` — different arity and
argument order. Not a Schwinger.jl bug, but if the two packages are meant to feel unified,
aligning the signature would help. (Out of scope for Schwinger.jl alone.)

---

## 4. Quick-win summary

| # | Item | Type | Effort | Payoff |
|---|------|------|--------|--------|
| 2.3 | Generic `energycurrents`/`chargecurrents`/`momentumdensities` for ED/ITensors | 🧩 | low | ED cross-checks |
| 2.9 | Exported `rehost`/`quench(state, Hnew)` | 🧩 | low | deletes `_withH` everywhere |
| 2.8 | State `normalize` + `savestate`/`loadstate` | 🧩 | low | removes `.psi` reach + reload boilerplate |
| 2.10 | Built-in observable menu for `evolve` | 🧩 | low | every driver stops re-writing it |
| 2.7 | `groupvelocity` / `dispersion` helper | 🧩 | low | every wavepacket script |
| 1.1 | Bond-consistent `energy_densities` (or loud docs) | 🐞 | low-med | correctness of energy continuity |
| 3.1 | Uniform `loweststates` return / accessors | ✂️ | low | removes a silent footgun |
| 2.1 | `act(local_op, ::WindowMPS)` on window tensors | 🧩 | med-high | removes the correlator tensor surgery |
| 2.2 | Exported real-time n-point correlator + E₀ phase | 🧩 | high | retires `correlators/correlators.jl` internals |

The three genuinely hard ones (2.1, 2.2, and to a lesser extent 2.5/2.6) all trace back to
one root cause: **local operators can't be applied to a `WindowMPS`.** Fixing 2.1 unlocks
most of the rest.

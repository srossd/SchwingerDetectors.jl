# SchwingerDetectors.jl

Real-time **detector correlators** for the Hamiltonian lattice Schwinger model, built as a
thin, documented layer over exported [`Schwinger.jl`](https://github.com/srossd/Schwinger.jl)
calls. A "detector" is a conserved current (`ChargeCurrent` j¹, `EnergyCurrent` 𝒥 = T⁰¹) or a
density, read at a fixed lattice site as a function of time.

Units are `g = 1` throughout: `ag = a·g` is the spacing, `mg = m/g` the mass.

## Installation

`SchwingerDetectors`, `Schwinger`, and Schwinger's dependency `MPSKitLEMPO` are all
unregistered, so add the three in a single `Pkg` operation:

```julia
using Pkg
Pkg.add([
    PackageSpec(url = "https://github.com/benjamints/MPSKitLEMPO.jl"),
    PackageSpec(url = "https://github.com/srossd/Schwinger.jl"),
    PackageSpec(url = "https://github.com/srossd/SchwingerDetectors.jl"),
])
```

## Usage

```julia
using Schwinger, SchwingerDetectors

# A finite θ = 0 lattice, then a sudden ±q probe-charge "string" quench.
m = build_model(; N = 64, ag = 0.2, mg = 0.0)
q = string_quench(m; length_links = 8)
c = 64 ÷ 2

# 1-point: watch the energy current 𝒥 radiate out to two fixed detectors.
t, series = detector_1pt(q.state, 6.0, [c - 8, c + 8];
                         kind = :energycurrent, nsteps = 60, maxbond = 128)

# 2-point: real-time ⟨j¹(x_D,t) j¹(x_S,0)⟩ on the vacuum (an eigenstate reference).
lat = m.lat
t2, C = detector_2pt(q.gs,
                     ChargeCurrent(lat, c - 8; backend = :MPSKit),   # source j¹ at t = 0
                     ChargeCurrent(lat, c + 8; backend = :MPSKit),   # detector j¹ at t
                     0.0:0.1:6.0; connected = true, maxlinkdim = 128)
```

More end-to-end scripts (each honours `SMOKE=1` for a seconds-long wiring run) live in
`demo.jl` and `examples/`.

## Capabilities

### State preparation
- `build_model` / `prepare_groundstate` — bundle a lattice + Hamiltonian; take the vacuum.
- `string_quench` — insert a static ±q probe-charge pair as a `θ/2π` step (keeps per-bond
  detectors working).
- `wilson_line_quench` — apply a gauge-invariant χ†…χ string to the vacuum (a meson at θ=0; a
  soliton–antisoliton bubble at θ=π).
- `local_operator_quench` — apply any exported local operator to a state.
- `moving_soliton` / `moving_meson` / `moving_wavepackets` — momentum-eigenstate
  quasiparticles Gaussian-enveloped into `WindowMPS` wavepackets (need `N = Inf`);
  `aligned_support` snaps a window to the unit cell.

### Evolve and measure
- `evolve_and_measure` — real-time-evolve while recording a menu of densities per TDVP step,
  with optional vacuum subtraction and checkpointing (safe on finite lattices and `WindowMPS`).
- `times` / `density_map` — read times and a `nt × nx` spacetime density array back out.
- `vacuum_tile` — tile an infinite-vacuum density to use as a subtraction background.

### Detectors
- `detector_1pt` — ⟨O(x_D, t)⟩ via the vectorized density; works on finite lattices **and**
  `WindowMPS`. `detector_1pt_operator` is the `expectation`-based variant (finite only).
- `detector_2pt` — real-time ⟨vac| O(x_D,t) O'(x_S,0) |vac⟩ on a finite lattice with an
  eigenstate reference (a named wrapper over `Schwinger.correlator2pt`).
- `detector_2pt_equal_time` — connected equal-time ⟨O O'⟩ *inside* a non-eigenstate quench
  state, evaluated on each snapshot.

## Caveats

1. Per-bond/site current operators (`ChargeCurrent`/`EnergyCurrent` via `act`/`expectation`)
   apply only to **plain finite-lattice** states. On a `WindowMPS` use the vectorized
   densities (`chargecurrents`, `energycurrents(…; pad=true)`, …); `detector_1pt` and
   `evolve_and_measure` do this and are safe on both, while the operator routines error
   clearly on a window state.
2. `detector_2pt` / `correlator2pt` subtract the vacuum phase e^{iE₀t} with the single number
   E₀ = ⟨vac|H|vac⟩ — exact only for an eigenstate reference. For a two-time correlator inside
   a quench state, use `detector_2pt_equal_time` on snapshots.

## Documentation

Full API docs (Documenter.jl) live in `docs/`; build locally with

```bash
julia --project=docs -e 'using Pkg; Pkg.develop([
    PackageSpec(url="https://github.com/benjamints/MPSKitLEMPO.jl"),
    PackageSpec(url="https://github.com/srossd/Schwinger.jl"),
    PackageSpec(path=".")]); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

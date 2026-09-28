# State preparation

Build a model, then pick a preparation. Every routine is exported `Schwinger.jl` API
(`Lattice`, `Hamiltonian`, `groundstate`, `loweststates`, `wavepacket`, `WilsonLine`,
`FermionField`, `act`, and the sudden-quench primitive `rehost`/`quench`).

```julia
m  = build_model(; N = 128, ag = 0.2, mg = 0.0)     # finite; N = Inf for wavepackets
gs = prepare_groundstate(m; energy_tol = 1e-9)
```

Preparations:

- **Static probe-charge string** — [`string_quench`](@ref) inserts a ±q pair as a `θ/2π`
  step (the Gauss-law field of a defect pair) on a *plain* N-site lattice, so per-bond
  detectors keep working. Uses the equivalent θ-step rather than an inert `DefectCharge`
  site, which per-bond operators cannot address.
- **Wilson-line quench** — [`wilson_line_quench`](@ref) applies the gauge-invariant
  ``\chi^\dagger\!\cdots\chi`` string to the vacuum (a meson at ``\theta=0``; a
  soliton–antisoliton bubble at ``\theta=\pi``, where `conjugate=true` picks the tensionless
  partner).
- **Generic local operator** — [`local_operator_quench`](@ref) for any exported local `op`.
- **Moving wavepackets** — [`moving_soliton`](@ref) / [`moving_meson`](@ref) build a
  momentum-eigenstate quasiparticle and Gaussian-envelope it into a `WindowMPS` (need
  `N = Inf`); [`moving_wavepackets`](@ref) places several disjoint packets (e.g. a colliding
  pair). [`aligned_support`](@ref) snaps a site window to the unit cell.

```@docs
build_model
prepare_groundstate
string_quench
wilson_line_quench
local_operator_quench
moving_soliton
moving_meson
moving_wavepackets
aligned_support
```

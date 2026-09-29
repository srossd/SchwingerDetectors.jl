# SchwingerDetectors.jl

Real-time **detector correlators** for the Hamiltonian lattice Schwinger model, built as
thin, documented glue over exported [`Schwinger.jl`](https://github.com/srossd/Schwinger.jl)
calls. A "detector" is a conserved current (`ChargeCurrent` ``j^1``, `EnergyCurrent`
``\mathcal{J}=T^{01}``) or a density, read at a fixed lattice site as a function of time.

Everything lives in the module `SchwingerDetectors`, spread across three files:

| file | concern | key functions |
|------|---------|---------------|
| `state_prep.jl` | prepare initial states | [`build_model`](@ref), [`prepare_groundstate`](@ref), [`string_quench`](@ref), [`wilson_line_quench`](@ref), [`local_operator_quench`](@ref), [`moving_soliton`](@ref), [`moving_meson`](@ref), [`moving_wavepackets`](@ref) |
| `measure.jl` | evolve + record densities | [`evolve_and_measure`](@ref), [`times`](@ref), [`density_map`](@ref), [`vacuum_tile`](@ref) |
| `detectors.jl` | detector 1-pt / 2-pt | [`detector_1pt`](@ref), [`detector_1pt_operator`](@ref), [`detector_2pt`](@ref), [`detector_2pt_equal_time`](@ref) |
| `charge_transfer.jl` | integrated-current correlator | [`charge_transfer_correlator`](@ref) |

## Table of contents

```@contents
Pages = ["man/stateprep.md", "man/measure.md", "man/detectors.md", "man/examples.md"]
Depth = 3
```

## Installation

`SchwingerDetectors` and its dependencies (`Schwinger`, and Schwinger's own dependency
`MPSKitLEMPO`) are unregistered, so add all three in one `Pkg` operation:

```julia
using Pkg
Pkg.add([
    PackageSpec(url = "https://github.com/benjamints/MPSKitLEMPO.jl"),
    PackageSpec(url = "https://github.com/srossd/Schwinger.jl"),
    PackageSpec(url = "https://github.com/srossd/SchwingerDetectors.jl"),
])
```

```julia
using SchwingerDetectors
```

Units are `g = 1` throughout, so `ag = a·g` is the spacing and `mg = m/g` the mass.

## When to use what

- **1-point** ``\langle O(x_D,t)\rangle`` — [`detector_1pt`](@ref) reads the vectorized
  density and works on both finite lattices and `WindowMPS` wavepackets;
  [`detector_1pt_operator`](@ref) is the `expectation`-based variant (finite only).
- **2-point, eigenstate reference** ``\langle\text{vac}|O(x_D,t)O'(x_S,0)|\text{vac}\rangle``
  — [`detector_2pt`](@ref), a named wrapper over `Schwinger.correlator2pt`.
- **2-point inside a quench (non-eigenstate) state** — [`detector_2pt_equal_time`](@ref) on
  each snapshot; the single-phase `correlator2pt` trick does not apply here.
- **Time-integrated current (charge-transfer) correlator in a non-eigenstate** —
  [`charge_transfer_correlator`](@ref) for ``\langle\psi|Q(b_L)Q(b_R)|\psi\rangle`` with
  ``Q(x)=\int_0^T j^1(x,t)dt``; the full two-time evolution is done explicitly.

## Caveats

1. Per-bond/site current operators (`ChargeCurrent`/`EnergyCurrent` via `act`/`expectation`)
   apply only to **plain finite-lattice** states. On a `WindowMPS` use the vectorized
   densities (`chargecurrents`, `energycurrents(…; pad=true)`, …); [`detector_1pt`](@ref) and
   [`evolve_and_measure`](@ref) do this and are safe on both. The operator routines error
   clearly on a window state.
2. `correlator2pt`/[`detector_2pt`](@ref) subtract the vacuum phase ``e^{iE_0 t}`` with the
   single number ``E_0=\langle\text{vac}|H|\text{vac}\rangle`` — exact only for an eigenstate
   reference. For a two-time correlator inside a quench state, use the equal-time route.

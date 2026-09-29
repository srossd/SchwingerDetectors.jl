# Detectors

A detector sits at a fixed lattice location and reads a conserved-current operator (or
density) versus time.

## One-point ``\langle O(x_D,t)\rangle``

[`detector_1pt`](@ref) evolves the state and returns the density `kind` at each fixed
location. It reads the vectorized density, so it is **universal** — finite lattices and
`WindowMPS` wavepackets alike:

```julia
t, series = detector_1pt(state, T, [xL, xR]; kind = :energycurrent, nsteps = 100,
                         two_site = true, maxbond = 128)
```

[`detector_1pt_operator`](@ref) is the `expectation`-based variant (finite lattices only),
matching the historical `step4c_detector.jl` measurement.

## Two-point

[`detector_2pt`](@ref) is a named wrapper over `Schwinger.correlator2pt`: the real-time
two-time correlator ``C(t)=\langle\text{vac}|\,O(x_D,t)\,O'(x_S,0)\,|\text{vac}\rangle`` on a
**finite lattice** with an **eigenstate** (vacuum) reference. The wrapper takes
`(source, detector)` order and returns `(tgrid, C)`:

```julia
lat = m.lat
t, C = detector_2pt(vac,
                    ChargeCurrent(lat, bS; backend = :MPSKit),   # source at t = 0
                    ChargeCurrent(lat, bD; backend = :MPSKit),   # detector at t
                    collect(0.0:0.1:10.0); connected = true, two_site = true, maxlinkdim = 128)
```

For the equal-time ``\langle j^1 j^1\rangle`` *inside* a non-eigenstate quench state as it
evolves, use [`detector_2pt_equal_time`](@ref) on each snapshot — the single-phase
`correlator2pt` subtraction does not apply to a non-eigenstate.

```@docs
detector_1pt
detector_1pt_operator
detector_2pt
detector_2pt_equal_time
```

## Time-integrated ("charge-transfer") correlator

For the correlator of the *time-integrated* current detectors ``Q(x)=\int_0^T j^1(x,t)\,dt``
in an arbitrary state (e.g. a Wilson-line quench, which is **not** an energy eigenstate),

```math
C(T) = \langle\psi|\,Q(b_L)\,Q(b_R)\,|\psi\rangle
     = \int_0^T\!\!\int_0^T dt_1\,dt_2\,\langle\psi|\,j^1(b_L,t_1)\,j^1(b_R,t_2)\,|\psi\rangle ,
```

use [`charge_transfer_correlator`](@ref). This is a genuinely two-time object: because
``|\psi\rangle`` is not stationary, the ``e^{\pm iHt}`` are real many-body evolutions and
**cannot** be replaced by a single scalar phase ``e^{iE_0 t}`` (the shortcut `correlator2pt`
uses, valid only for an eigenstate reference). Since ``Q`` is Hermitian the correlator equals
``\langle u|v\rangle`` with ``|v\rangle=Q(b_R)|\psi\rangle``, ``|u\rangle=Q(b_L)|\psi\rangle``;
the overall ``e^{iHT}`` cancels, so the result is phase-convention-free.

Two independent algorithms are provided and cross-checked (`method = :both`): an O(nsteps)
forward accumulation (`:accumulated`) and an O(nsteps²) reference built from the full two-time
matrix (`:matrix`). The return value is the **cumulative curve** ``C(t_m)`` at every grid
time, i.e. a built-in convergence scan in the upper limit ``T``, with both `full` and
`connected` variants.

```@docs
charge_transfer_correlator
```

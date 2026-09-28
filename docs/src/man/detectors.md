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

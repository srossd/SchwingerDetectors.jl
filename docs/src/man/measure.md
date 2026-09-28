# Evolve and measure

[`evolve_and_measure`](@ref) is a thin convenience over `Schwinger.evolve` +
`standard_densities`: real-time-evolve a state and record a menu of densities at every TDVP
step, with optional vacuum subtraction and checkpointing. All densities are exported
observables and work on finite-lattice states **and** `WindowMPS` wavepackets.

```julia
evolved, obs = evolve_and_measure(state, T; nsteps = 100, two_site = true, maxbond = 128,
                                  densities = [:energy, :electricfield, :current,
                                               :energycurrent, :pseudoscalar])
t   = times(obs)                        # Vector of times
Jxt = density_map(obs, :energycurrent)  # nt × nx spacetime array
```

The density menu (`Schwinger.standard_densities`): `:charge` (``\rho=j^0``), `:occupation`,
`:electricfield`, `:energy`, `:current` (``j^1`` per bond), `:energycurrent`
(``\mathcal{J}`` per site), `:pseudoscalar`, `:momentum`. For the ``\mathcal{J}``-consistent
energy density call `energy_densities(state; convention=:bond)` directly.

- `subtract = Dict(:energy => vacuum_tile(vac, :energy, W), …)` stores only the excess over
  an infinite-lattice vacuum; [`vacuum_tile`](@ref) tiles the staggered background.
- `save`/`save_every` checkpoint; `grow = true` enables adaptive `WindowMPS` growth.

```@docs
evolve_and_measure
times
density_map
vacuum_tile
```

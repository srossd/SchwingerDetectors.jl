# =============================================================================
# measure.jl — time-evolve while measuring densities.
#
# Thin conveniences over the built-in `evolve` + `standard_densities` observable menu
# (Schwinger.jl now ships the menu itself). `evolve_and_measure` adds vacuum-subtraction,
# checkpoint-saving, and bond-dimension ergonomics on top; `density_map`/`times`/
# `vacuum_tile` help read the result back. All density functions are exported observables
# and work on plain finite-lattice states AND WindowMPS wavepackets.
# =============================================================================

# Built-in menu names (from Schwinger.standard_densities): :charge, :occupation,
# :electricfield, :energy, :current (j¹), :energycurrent (𝒥), :pseudoscalar, :momentum.

"""
    evolve_and_measure(state, T; nsteps, densities=[:energy,:electricfield,:current,:pseudoscalar],
                       two_site=true, maxbond=nothing, subtract=nothing,
                       save=nothing, save_every=1, grow=nothing, evolve_kwargs=(;))
        -> (evolved_state, obs)

Real-time-evolve `state` to time `T` in `nsteps` TDVP steps, recording each density in
`densities` (keys of the built-in `standard_densities` menu) at every step. Returns the
evolved state and the `Observers` object; read a field with [`times`] / [`density_map`].

- `two_site` + `maxbond`: two-site TDVP capped at bond dim `maxbond`; `two_site=false`
  uses single-site TDVP (fixed bond dim).
- `subtract`: optional `Dict(name => vacuum_vector)` subtracted from that density each
  step (e.g. tile the infinite-vacuum background so only the excess is stored — see
  [`vacuum_tile`]).
- `save`: optional `(state, obs) -> nothing` checkpoint, called every `save_every` steps.
- `grow=true`: adaptive window growth for a WindowMPS (forwarded to `evolve`).
"""
function evolve_and_measure(state, T::Real; nsteps::Int,
                            densities = [:energy, :electricfield, :current, :pseudoscalar],
                            two_site::Bool = true, maxbond = nothing,
                            subtract = nothing, save = nothing, save_every::Int = 1,
                            grow = nothing, evolve_kwargs = (;))
    base = standard_densities(densities)                 # built-in menu: name => (ψ,t)->vector
    obsdict = Dict{String,Function}()
    for (name, f) in base
        sub = isnothing(subtract) ? nothing : get(subtract, Symbol(name), nothing)
        obsdict[name] = isnothing(sub) ? f : ((ψ, t) -> f(ψ, t) .- sub)
    end
    ckpt = isnothing(save) ? nothing : ((ψ, t, step, obs) -> save(ψ, obs))
    mb   = two_site ? maxbond : nothing                  # single-site TDVP ignores maxlinkdim
    return evolve(state, T; nsteps = nsteps, two_site = two_site, maxlinkdim = mb,
                  observable = obsdict, checkpoint = ckpt, checkpoint_every = save_every,
                  grow = grow, evolve_kwargs...)
end

# -----------------------------------------------------------------------------
# convenience readers
# -----------------------------------------------------------------------------

"Vector of recorded times."
times(obs) = collect(obs.time)

"""
    density_map(obs, name) -> Matrix (nt × nx)

Stack a recorded density into a spacetime array: row = timestep, column = site/bond.
`name` is the density symbol/string used in `densities`. Assumes F=1 (vector per step).
"""
function density_map(obs, name)
    rows = collect(getproperty(obs, Symbol(name)))     # Vector (over t) of Vectors (over x)
    return permutedims(reduce(hcat, rows))             # nt × nx
end

"""
    vacuum_tile(vac, density, nsites) -> Vector

Tile a density measured on an INFINITE-lattice vacuum `vac` (length = unit cell) out to
`nsites`, respecting the staggered sublattice. Use as a `subtract` entry so a wavepacket
run stores only the excess over vacuum. `density` is a built-in menu symbol.
"""
function vacuum_tile(vac, density::Symbol, nsites::Int)
    f = standard_densities([density])[String(density)]
    base = f(vac, 0.0)
    L = length(base)
    return [base[mod1(i, L)] for i in 1:nsites]
end

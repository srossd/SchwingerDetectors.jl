# Examples

Three runnable scripts sit at the repo root / `examples/`. Each honours `SMOKE=1` for a
tiny, seconds-long wiring run and reads its parameters from environment variables.

## `demo.jl` — end-to-end smoke

Builds a string quench, reads an energy-current 1-pt series while recording density maps, and
computes a vacuum ``\langle j^1 j^1\rangle_c`` 2-pt correlator — exercising all three files.

```bash
SMOKE=1 julia --project=. demo.jl
AG=0.2 MG=0.0 N=64 julia --project=. demo.jl
```

## `examples/energy_1pt_probe_charges.jl` — energy 1-pt

Insert a static ±q probe-charge pair (a θ-step [`string_quench`](@ref)) and watch the energy
current ``\mathcal{J}=T^{01}`` radiate outward at two fixed detectors — the one-point function
``\langle\mathcal{J}(x_D,t)\rangle``. Also stores the bond-consistent energy density
(`energy_densities(…; convention=:bond)`), the density that partners ``\mathcal{J}`` in
``\partial_t h_n = \mathcal{J}_n - \mathcal{J}_{n+1}``.

```bash
SMOKE=1 julia --project=. examples/energy_1pt_probe_charges.jl
AG=0.2 MG=0.0 N=128 STRING_L=16 DET=24 T=14 julia --project=. examples/energy_1pt_probe_charges.jl
```

## `examples/charge_charge_wilson.jl` — charge–charge 2-pt inside a quench

Create a charge pair joined by a flux string, ``|\psi\rangle = W|\text{vac}\rangle``
([`wilson_line_quench`](@ref)), and measure the connected equal-time correlator
``C(R,t)=\langle\psi(t)|\,j^1(c+R)\,j^1(c-R)\,|\psi(t)\rangle_c`` on every TDVP snapshot with
[`detector_2pt_equal_time`](@ref). Because ``|\psi\rangle`` is not an eigenstate, this
in-state expectation — not `correlator2pt` — is the correct object.

```bash
SMOKE=1 julia --project=. examples/charge_charge_wilson.jl
AG=0.2 MG=0.0 N=128 NL=16 T=12 RLIST=8,16,24 julia --project=. examples/charge_charge_wilson.jl
```

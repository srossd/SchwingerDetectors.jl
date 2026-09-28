# =============================================================================
# SchwingerToolkit — consolidated real-time Schwinger-model capabilities.
#
# A small library over EXPORTED Schwinger.jl calls, grouped into three concerns:
#   • state_prep.jl — prepare quench / wavepacket initial states
#   • measure.jl    — time-evolve while recording densities
#   • detectors.jl  — detector-operator 1-point and 2-point functions
#
# Usage (from the project root, environment active):
#   include("Schwinger/toolkit/SchwingerToolkit.jl")
#   using .SchwingerToolkit
# See README.md for the capability map and the exported-API caveats.
# =============================================================================
module SchwingerToolkit

using Schwinger
using MPSKit
using LinearAlgebra

include("state_prep.jl")
include("measure.jl")
include("detectors.jl")

# --- model / state preparation (rehost/quench, vacuumof, quasiparticle come from Schwinger) ---
export build_model, prepare_groundstate
export string_quench, wilson_line_quench, local_operator_quench
export moving_soliton, moving_meson, moving_wavepackets, aligned_support

# --- evolve + measure (standard_densities comes from Schwinger) ---
export evolve_and_measure, times, density_map, vacuum_tile

# --- detectors (correlator2pt comes from Schwinger; detector_2pt wraps it) ---
export detector_1pt, detector_1pt_operator, detector_2pt, detector_2pt_equal_time

end # module

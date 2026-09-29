# Shared configuration for the charge-transfer production sweep.
# All lengths in units g = 1 (ag = a·g the spacing).
using Schwinger, SchwingerDetectors

# Defaults are the ag=0.2, same-physical-box (L=102.4/g) plan; override any via ENV.
const N       = parse(Int,     get(ENV, "N",        "512"))    # finite lattice sites
const AG      = parse(Float64, get(ENV, "AG",       "0.2"))    # a·g  (physical box L = N·ag = 102.4)
const MG      = parse(Float64, get(ENV, "MG",       "1.0"))    # m/g
const GX_WL   = parse(Float64, get(ENV, "GX_WL",    "2.0"))    # Wilson-line length in g·x  -> 10 sites
const RLIST   = [parse(Int, x) for x in split(get(ENV, "RLIST", "10,15,20"), ",")]  # detector offsets (sites) = phys 2,3,4
const T       = parse(Float64, get(ENV, "T",        "90.0"))   # evolution time (g·t); boundary-clean to ~98
const DT      = parse(Float64, get(ENV, "DT",       "0.2"))    # grid spacing (g·t); benchmarked
const SUBSTEP = parse(Int,     get(ENV, "SUBSTEPS", "1"))      # TDVP substeps per grid gap
const MAXBOND = parse(Int,     get(ENV, "MAXBOND",  "256"))    # evolution bond-dim cap
const GS_BOND = parse(Int,     get(ENV, "GS_BOND",  "96"))     # groundstate bond-dim (gapped -> small)
const GS_TOL  = parse(Float64, get(ENV, "GS_TOL",   "1e-6"))   # groundstate energy tol
const HEAT_STRIDE = parse(Int, get(ENV, "HEAT_STRIDE", "5"))   # record density maps every k steps
const EAVG_TOL = parse(Float64, get(ENV, "EAVG_TOL", "0.5"))   # |<E_avg>| must be below this (+eps)

# θ/π  ->  θ2π, folded into (-0.5, 0.5].  θ is periodic mod 2π (θ2π mod 1), so e.g. θ/π=1.2
# (θ2π=0.6) is the SAME physics as θ2π=-0.4.  The finite DMRG starts from a random L≈0 state
# with no branch bias, so it lands on the TRUE (screened, |<E>|<0.5) vacuum only when θ2π is in
# (-0.5,0.5] — confirmed by energy (θ2π=-0.4: E0=-308.22 < false +0.6 branch E0=-307.43).
theta2pi(theta_over_pi) = (r = theta_over_pi / 2; r - round(r))

# center bond/site of the lattice and the Wilson-line endpoints (gx = GX_WL wide, centered)
const CENTER = N ÷ 2
wilson_endpoints() = (CENTER - round(Int, GX_WL / AG) ÷ 2, CENTER - round(Int, GX_WL / AG) ÷ 2 + round(Int, GX_WL / AG))

model(theta2pi_val) = build_model(; N = N, F = 1, q = 1, ag = AG, mg = MG, theta2pi = theta2pi_val)

# a filesystem-safe tag for a θ/π value, e.g. 1.2 -> "1p2"
thetatag(theta_over_pi) = replace(string(theta_over_pi), "." => "p")

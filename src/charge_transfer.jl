# =============================================================================
# charge_transfer.jl — correlator of time-integrated current ("charge-transfer")
# detectors in an arbitrary (non-eigenstate) state.
#
# The detector is the net charge transported across a bond over [0, T]:
#
#     Q(x) = ∫₀ᵀ j¹(x, t) dt ,   j¹(x, t) = e^{iHt} j¹(x) e^{-iHt}   (Heisenberg).
#
# We compute the two-detector correlator in a state |ψ⟩ (e.g. a Wilson-line quench,
# which is NOT an energy eigenstate):
#
#     C(T) = ⟨ψ| Q(x_L) Q(x_R) |ψ⟩
#          = ∫₀ᵀ dt₁ ∫₀ᵀ dt₂ ⟨ψ| j¹(x_L, t₁) j¹(x_R, t₂) |ψ⟩ .
#
# Because |ψ⟩ is not stationary, the two time evolutions e^{±iHt} are GENUINE
# many-body evolutions — they can NOT be collapsed to a single scalar phase e^{iE₀t}
# (that trick, used by `correlator2pt`, is exact only for an eigenstate reference).
# The `E₀` phase is therefore irrelevant here: the overall e^{iHT} cancels between the
# bra and ket (see below), so the result is convention-independent and nothing is
# "subtracted". Only the physical relative phases between energy components survive,
# and they are captured by doing the actual evolution.
#
# Since Q(x) is Hermitian, C(T) = ⟨u|v⟩ with |v⟩ = Q(x_R)|ψ⟩, |u⟩ = Q(x_L)|ψ⟩, and
#     |v⟩ = ∫₀ᵀ e^{iHt} j¹(x_R) e^{-iHt} |ψ⟩ dt .
#
# Two independent algorithms are implemented and cross-checked:
#   :accumulated  — O(M) forward pass. Build |u⟩,|v⟩ incrementally with a running
#                   (Riemann) accumulator plus trapezoid end corrections; C(T) = ⟨u|v⟩.
#   :matrix       — O(M²) reference. Backward-evolve each kicked snapshot to a common
#                   reference time to form the full 2-time matrix D(t₁,t₂), then
#                   2-D trapezoid. Transparent but scales poorly; use for validation.
#
# Both return the whole cumulative curve C(t_m) for every grid time t_m ≤ T (a free
# convergence scan in the upper limit T). Finite-lattice states only (per-bond `act`).
# =============================================================================

# add SchwingerStates by summing the underlying MPS (bond dims add; a subsequent
# `evolve` retruncates). Used only for the accumulators.
_addstates(ss...) = MPSKitState(ss[1].hamiltonian, sum(s.psi for s in ss))

# cumulative composite-trapezoid weights on the grid 0 = t₀ < … < t_m (uniform Δt).
function _trapweights(m::Int, dt::Float64)
    m == 0 && return [0.0]
    w = fill(dt, m + 1); w[1] = dt / 2; w[end] = dt / 2
    return w
end

# forward pass: snapshots |ψ(t_k)⟩ and the one-point currents ⟨j_bL⟩(t_k), ⟨j_bR⟩(t_k).
function _forward_snapshots(state, jL, jR, M, dt, substeps, two_site, maxbond)
    mb = two_site ? maxbond : nothing
    snaps = Vector{MPSKitState}(undef, M + 1); snaps[1] = state
    jLt = zeros(Float64, M + 1); jRt = zeros(Float64, M + 1)
    jLt[1] = real(expectation(jL, state)); jRt[1] = real(expectation(jR, state))
    for k in 1:M
        snaps[k + 1] = evolve(snaps[k], dt; nsteps = substeps, two_site = two_site,
                              maxlinkdim = mb)[1]
        jLt[k + 1] = real(expectation(jL, snaps[k + 1]))
        jRt[k + 1] = real(expectation(jR, snaps[k + 1]))
    end
    return snaps, jLt, jRt
end

# cumulative time integrals of the one-point currents, ∫₀^{t_m} ⟨j⟩ dt.
function _cumintegral(f, M, dt)
    out = zeros(Float64, M + 1)
    for m in 1:M
        w = _trapweights(m, dt); out[m + 1] = sum(w .* f[1:m + 1])
    end
    return out
end

# ---- :matrix (O(M²) reference) ------------------------------------------------
function _qq_matrix(snaps, jL, jR, M, dt, substeps, two_site, maxbond)
    mb = two_site ? maxbond : nothing
    bwd(s, t) = evolve(s, -t; nsteps = max(1, substeps), two_site = two_site, maxlinkdim = mb)[1]
    ts = collect(0:M) .* dt
    A = Vector{MPSKitState}(undef, M + 1); B = Vector{MPSKitState}(undef, M + 1)
    for k in 0:M
        ak = act(jR, snaps[k + 1]); bk = act(jL, snaps[k + 1])
        A[k + 1] = k == 0 ? ak : bwd(ak, ts[k + 1])   # e^{iH t_k} jR |ψ(t_k)⟩
        B[k + 1] = k == 0 ? bk : bwd(bk, ts[k + 1])   # e^{iH t_k} jL |ψ(t_k)⟩
    end
    D = [dot(B[j + 1], A[k + 1]) for j in 0:M, k in 0:M]   # ⟨ψ| jL(t_j) jR(t_k) |ψ⟩
    full = zeros(ComplexF64, M + 1)
    for m in 0:M
        w = _trapweights(m, dt); s = 0.0 + 0im
        for j in 0:m, k in 0:m; s += w[j + 1] * w[k + 1] * D[j + 1, k + 1]; end
        full[m + 1] = s
    end
    return full
end

# ---- :accumulated (O(M) forward) ---------------------------------------------
function _qq_accumulated(snaps, jL, jR, M, dt, substeps, two_site, maxbond)
    mb = two_site ? maxbond : nothing
    fwd(s) = evolve(s, dt; nsteps = max(1, substeps), two_site = two_site, maxlinkdim = mb)[1]
    a0 = act(jR, snaps[1]); b0 = act(jL, snaps[1])
    Rv = dt * a0; Ru = dt * b0            # Riemann accumulators (weight Δt)
    e0v = a0;     e0u = b0                 # e^{-iH t_k} α₀ / β₀ (for endpoint correction)
    full = zeros(ComplexF64, M + 1)        # full[1] = C(0) = 0
    for k in 1:M
        ak = act(jR, snaps[k + 1]); bk = act(jL, snaps[k + 1])
        Rv = _addstates(fwd(Rv), dt * ak); e0v = fwd(e0v)
        Ru = _addstates(fwd(Ru), dt * bk); e0u = fwd(e0u)
        # trapezoid = Riemann − ½Δt(endpoints); endpoints are α_k and e^{-iH t_k} α₀
        vk = _addstates(Rv, (-dt / 2) * ak, (-dt / 2) * e0v)
        uk = _addstates(Ru, (-dt / 2) * bk, (-dt / 2) * e0u)
        full[k + 1] = dot(uk, vk)          # e^{iH t_k} cancels between ⟨u| and |v⟩
    end
    return full
end

"""
    charge_transfer_correlator(state, T, (bL, bR); nsteps, substeps=1, two_site=true,
                               maxbond=nothing, method=:both, cross_tol=1e-2)
        -> NamedTuple

Correlator of the time-integrated current ("charge-transfer") detectors
`Q(x) = ∫₀ᵀ j¹(x,t) dt` on a FINITE lattice, in an arbitrary state `state` (typically a
non-eigenstate such as a Wilson-line quench):

    C(T) = ⟨state| Q(bL) Q(bR) |state⟩
         = ∫₀ᵀ dt₁ ∫₀ᵀ dt₂ ⟨state| j¹(bL,t₁) j¹(bR,t₂) |state⟩ .

`bL, bR` are bond indices (e.g. `c-R, c+R` about a centre `c`). `state` should be
normalized (`Q` is linear in it). The evolution grid is `t_k = k·T/nsteps`, `k=0…nsteps`;
`substeps` TDVP substeps are taken per gap; two-site TDVP is capped at `maxbond`.

Returns `(; t, full, connected, QL, QR, method, discrepancy, full_matrix, full_accumulated)`
where each of `full`/`connected` is the CUMULATIVE curve `C(t_m)` at every grid time
`t_m ≤ T` — i.e. a built-in convergence scan in the upper limit. `connected` subtracts the
disconnected piece `⟨Q(bL)⟩⟨Q(bR)⟩`; `QL`,`QR` are the cumulative `∫₀^{t_m}⟨j¹⟩dt`.

`method`:
- `:accumulated` — O(nsteps) forward pass (recommended for large `T`).
- `:matrix` — O(nsteps²) reference via the full two-time matrix (transparent, slow).
- `:both` (default) — run both and report `discrepancy = max|C_acc − C_mat|`; `full` is the
  accumulated result. A `discrepancy` above `cross_tol` (relative to `max|full|`) is a
  numerical-convergence warning — increase `maxbond`/`substeps`, do not ignore it.

Why no `E₀` phase: |state⟩ is not an eigenstate, so ⟨state|e^{iHt} is not a scalar times
⟨state|; the genuine two-time evolution is done explicitly. The overall e^{iHT} cancels in
`⟨Q(bL)state | Q(bR)state⟩`, so the result is phase-convention-independent.
"""
function charge_transfer_correlator(state::MPSKitState, T::Real, bonds::Tuple{Int,Int};
                                    nsteps::Int, substeps::Int = 1, two_site::Bool = true,
                                    maxbond = nothing, method::Symbol = :both,
                                    cross_tol::Real = 1e-2)
    _assert_not_window(state, "charge_transfer_correlator")
    method in (:both, :accumulated, :matrix) ||
        throw(ArgumentError("method must be :both, :accumulated, or :matrix"))
    lat = state.hamiltonian.lattice
    bL, bR = bonds
    jL = ChargeCurrent(lat, bL; backend = :MPSKit)
    jR = ChargeCurrent(lat, bR; backend = :MPSKit)

    M = nsteps; dt = float(T) / M
    ts = collect(0:M) .* dt
    snaps, jLt, jRt = _forward_snapshots(state, jL, jR, M, dt, substeps, two_site, maxbond)
    QL = _cumintegral(jLt, M, dt); QR = _cumintegral(jRt, M, dt)

    facc = (method === :accumulated || method === :both) ?
           _qq_accumulated(snaps, jL, jR, M, dt, substeps, two_site, maxbond) : nothing
    fmat = (method === :matrix || method === :both) ?
           _qq_matrix(snaps, jL, jR, M, dt, substeps, two_site, maxbond) : nothing

    full = method === :matrix ? fmat : facc
    disc = (method === :both) ? maximum(abs.(facc .- fmat)) : NaN
    if method === :both
        scale = max(maximum(abs.(full)), eps())
        disc / scale > cross_tol && @warn "charge_transfer_correlator: methods disagree by " *
            "$(round(disc/scale; sigdigits=3)) (rel) > cross_tol=$cross_tol; " *
            "increase maxbond/substeps." discrepancy = disc
    end
    connected = full .- (QL .* QR)
    return (; t = ts, full = full, connected = connected, QL = QL, QR = QR,
            method = method, discrepancy = disc,
            full_matrix = fmat, full_accumulated = facc)
end

"""
    charge_transfer_correlator(state, T, pairs::AbstractVector{<:Tuple{Integer,Integer}};
                               nsteps, substeps=1, two_site=true, maxbond=nothing,
                               observe=nothing) -> NamedTuple

Streaming, memory-safe, multi-detector version of the `:accumulated` method: computes the
charge-transfer correlator `C(T)=⟨state|Q(bL)Q(bR)|state⟩` for SEVERAL bond pairs at once in a
SINGLE forward evolution of `state`, WITHOUT storing snapshots (suitable for large lattices /
long times where the O(nsteps)-snapshot storage of the single-pair method is prohibitive).

`pairs` is a vector of `(bL, bR)` bond tuples. Each pair carries its own O(1) set of forward
accumulators (see the file header); the state itself is evolved once. Returns
`(; t, pairs, full, connected, QL, QR)` where `full`, `connected`, `QL`, `QR` are `Vector`s
indexed like `pairs`, each the cumulative curve over the grid `t`. `observe(ψ, t, k)` is an
optional callback run after each step on the evolved physical state (e.g. to record density
spacetime maps for heatmaps in the same pass). This is the routine used for production runs; it
reproduces the validated single-pair `:accumulated` result (they share the same recursion).
"""
function charge_transfer_correlator(state::MPSKitState, T::Real,
                                    pairs::AbstractVector{<:Tuple{Integer,Integer}};
                                    nsteps::Int, substeps::Int = 1, two_site::Bool = true,
                                    maxbond = nothing, observe = nothing)
    _assert_not_window(state, "charge_transfer_correlator")
    isempty(pairs) && throw(ArgumentError("pairs must be non-empty"))
    lat = state.hamiltonian.lattice
    P = length(pairs)
    jLs = [ChargeCurrent(lat, p[1]; backend = :MPSKit) for p in pairs]
    jRs = [ChargeCurrent(lat, p[2]; backend = :MPSKit) for p in pairs]

    M = nsteps; dt = float(T) / M; mb = two_site ? maxbond : nothing
    fwd(s) = evolve(s, dt; nsteps = substeps, two_site = two_site, maxlinkdim = mb)[1]

    ψ = state
    Rv = Vector{MPSKitState}(undef, P); Ru = similar(Rv)
    e0v = Vector{MPSKitState}(undef, P); e0u = similar(e0v)
    full = [zeros(ComplexF64, M + 1) for _ in 1:P]
    jLt = [zeros(Float64, M + 1) for _ in 1:P]; jRt = [zeros(Float64, M + 1) for _ in 1:P]
    for p in 1:P
        a0 = act(jRs[p], ψ); b0 = act(jLs[p], ψ)
        Rv[p] = dt * a0; Ru[p] = dt * b0; e0v[p] = a0; e0u[p] = b0
        jLt[p][1] = real(expectation(jLs[p], ψ)); jRt[p][1] = real(expectation(jRs[p], ψ))
    end
    observe !== nothing && observe(ψ, 0.0, 0)

    for k in 1:M
        ψ = fwd(ψ)
        for p in 1:P
            ak = act(jRs[p], ψ); bk = act(jLs[p], ψ)
            Rv[p] = _addstates(fwd(Rv[p]), dt * ak); e0v[p] = fwd(e0v[p])
            Ru[p] = _addstates(fwd(Ru[p]), dt * bk); e0u[p] = fwd(e0u[p])
            vk = _addstates(Rv[p], (-dt / 2) * ak, (-dt / 2) * e0v[p])
            uk = _addstates(Ru[p], (-dt / 2) * bk, (-dt / 2) * e0u[p])
            full[p][k + 1] = dot(uk, vk)
            jLt[p][k + 1] = real(expectation(jLs[p], ψ))
            jRt[p][k + 1] = real(expectation(jRs[p], ψ))
        end
        observe !== nothing && observe(ψ, k * dt, k)
    end

    ts = collect(0:M) .* dt
    QL = [_cumintegral(jLt[p], M, dt) for p in 1:P]
    QR = [_cumintegral(jRt[p], M, dt) for p in 1:P]
    connected = [full[p] .- (QL[p] .* QR[p]) for p in 1:P]
    return (; t = ts, pairs = collect(pairs), full = full, connected = connected,
            QL = QL, QR = QR)
end

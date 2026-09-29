# MPS persistence via Julia's Serialization stdlib. (Schwinger's savestate/loadstate are a newer
# feature not present in the GitHub-main build used on the cluster, so we roll our own: serialize
# the underlying MPS tensors and rewrap with a freshly built Hamiltonian on load. Same Julia +
# package versions on save/load, so this round-trips exactly.) Atomic write via tmp + rename.
using Serialization
using Schwinger: MPSKitState

save_state(path, st) = (tmp = path * ".tmp"; open(io -> serialize(io, st.psi), tmp, "w"); mv(tmp, path; force = true))
load_state(path, H)  = MPSKitState(H, open(deserialize, path))

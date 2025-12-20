module DFTatom

using GaussianBasis
using LinearAlgebra
using OMEinsum
using Printf

include("AufbauSelection.jl")
export get_basis_angular_momentum, select_orbitals_aufbau
include("HF.jl")
export init_conf, SCF
include("LDA.jl")
export KS_SCF, run_lda_c, run_lda_he

end

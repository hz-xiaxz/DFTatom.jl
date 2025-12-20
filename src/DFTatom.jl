module DFTatom

using GaussianBasis
using LinearAlgebra
using OMEinsum
using Printf

include("AufbauSelection.jl")
include("HF.jl")
export init_conf, SCF
include("LDA.jl")
export KS_SCF, run_lda_c, run_lda_he

end

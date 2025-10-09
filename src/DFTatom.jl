module DFTatom

using GaussianBasis
using LinearAlgebra
using OMEinsum
using Printf

include("HF.jl")
export init_conf, SCF
include("LDA.jl")
export run_lda_c, run_lda_he

end

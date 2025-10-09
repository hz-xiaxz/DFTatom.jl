import Pkg
Pkg.activate(".")

using DFTatom
using GaussianBasis
using Printf

# This script runs a Hartree-Fock calculation for a Carbon atom
# with different basis sets for the singlet state.

println("Calculating Carbon atom (singlet state, N_up=3, N_down=3)")

# List of basis sets to use, based on README.md
basis_sets = [
    "sto-3g",
    "sto-6g",
    "6-31g",
    "6-31g_st_", # Corresponds to 6-31g*
    "cc-pvdz",
]

# For this calculation, we use N_up=3, N_down=3, which corresponds to an
# excited singlet state for the Carbon atom (1s²2s²2p²).
N_up = 3
N_down = 3

println("| Basis Set | HF Energy (Hartree) |")
println("|:----------|:--------------------|")

for basis_name in basis_sets
    # Construct the basis set for a Carbon atom at the origin
    bset = BasisSet(basis_name, "C 0.0 0.0 0.0")

    # Increase maxiter for basis sets known to converge slowly
    max_iterations = 100
    if basis_name == "cc-pvdz"
        max_iterations = 3000
    end

    # Perform the SCF calculation
    try
        result = SCF(bset; N_up = N_up, N_down = N_down, maxiter = max_iterations, α=0.5)
        @printf "| %-9s | %-19.10f |\n" basis_name result.energy
    catch e
        @printf "| %-9s | %-19s |\n" basis_name "Calculation Failed"
        println(e)
    end
end

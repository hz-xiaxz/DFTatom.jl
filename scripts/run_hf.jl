import Pkg
Pkg.activate(".")

using DFTatom
using GaussianBasis
using Printf

# This script runs Hartree-Fock calculations for a Carbon atom
# with different basis sets for the GROUND STATE (triplet, ³P).

println("="^70)
println("Carbon Atom Ground State (³P, N_up=4, N_down=2)")
println("Unrestricted Hartree-Fock (UHF) with different basis sets")
println("="^70)

# List of basis sets to test
basis_sets = [
    "sto-3g",
    "sto-6g",
    "6-31g",
    "6-31g_st_", # Corresponds to 6-31g*
    "cc-pvdz",
]

# Ground state configuration for Carbon
N_up = 4
N_down = 2

println("\n| Basis Set | HF Energy (Hartree) | Status |")
println("|:----------|:--------------------|:-------|")

for basis_name in basis_sets
    local bset = BasisSet(basis_name, "C 0.0 0.0 0.0")

    try
        # Use Aufbau with favor_high_m=false to prefer pz for triplet
        result = SCF(
            bset;
            N_up=N_up,
            N_down=N_down,
            maxiter=200,
            α=0.7,
            use_aufbau=true,
            favor_high_m=false
        )
        @printf "| %-9s | %19.10f | %-6s |\n" basis_name result.energy "OK"
    catch e
        @printf "| %-9s | %-19s | %-6s |\n" basis_name "---" "FAILED"
        println("  Error: ", e)
    end
end

println("\n" * "="^70)
println("Note: This is the Carbon atom GROUND STATE (³P triplet)")
println("Configuration: 1s² 2s² 2p² with N_up=4, N_down=2")
println("="^70)

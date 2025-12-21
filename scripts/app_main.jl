# Entry point for standalone app
using DFTatom
using GaussianBasis
using Printf

# Include orbital analysis utilities
include(joinpath(@__DIR__, "OrbitalAnalysis.jl"))

function julia_main()::Cint
    try
        println("="^80)
        println("     Atomic Electronic Structure Calculations")
        println("     Hartree-Fock (UHF) and DFT-LDA (LSDA) Methods")
        println("="^80)

        # ============================================================================
        # HYDROGEN ATOM (Ground State: 1s¹, ²S)
        # ============================================================================

        println("\n\n")
        println("#" * "="^78 * "#")
        println("#  HYDROGEN ATOM (H)  ")
        println("#  Ground State: 1s¹ (²S)")
        println("#" * "="^78 * "#")

        # --- Hartree-Fock ---
        println("\n" * "-"^80)
        println("  Method: Unrestricted Hartree-Fock (UHF)")
        println("-"^80)

        bset_H = BasisSet("6-31g_st_", "H 0.0 0.0 0.0")
        result_H_HF = DFTatom.SCF(bset_H; N_up=1, N_down=0, maxiter=100, α=0.8)

        print_orbital_analysis(bset_H, result_H_HF, "UHF"; N_up=1, N_down=0)

        # Save results
        save_results("hydrogen_HF.txt", bset_H, result_H_HF, "UHF"; N_up=1, N_down=0)

        println("\n" * "-"^80)
        println("Note: LSDA is not applicable for single-electron systems like Hydrogen.")
        println("LSDA requires multiple electrons to properly describe exchange-correlation.")
        println("For H, only Hartree-Fock results are physically meaningful.")
        println("-"^80)


        # ============================================================================
        # CARBON ATOM (Ground State: 1s² 2s² 2p², ³P)
        # ============================================================================

        println("\n\n")
        println("#" * "="^78 * "#")
        println("#  CARBON ATOM (C)  ")
        println("#  Ground State: 1s² 2s² 2p² (³P)")
        println("#  Configuration: N_up=4, N_down=2 (triplet)")
        println("#" * "="^78 * "#")

        # --- Hartree-Fock ---
        println("\n" * "-"^80)
        println("  Method: Unrestricted Hartree-Fock (UHF)")
        println("-"^80)

        bset_C = BasisSet("6-31g_st_", "C 0.0 0.0 0.0")
        result_C_HF = DFTatom.SCF(
            bset_C;
            N_up=4,
            N_down=2,
            maxiter=100,
            α=0.8,
            use_aufbau=true,
            favor_high_m=false
        )

        print_orbital_analysis(bset_C, result_C_HF, "UHF"; N_up=4, N_down=2)

        # --- DFT-LDA ---
        println("\n" * "-"^80)
        println("  Method: DFT with Local Spin Density Approximation (LSDA)")
        println("-"^80)

        result_C_LDA = DFTatom.KS_SCF(
            bset_C;
            N_up=4,
            N_down=2,
            maxiter=1000,
            α=0.5,
            use_aufbau=true,
            favor_high_m=false,
            tol=1e-5
        )

        print_orbital_analysis(bset_C, result_C_LDA, "LSDA"; N_up=4, N_down=2)

        # Save results
        save_results("carbon_HF.txt", bset_C, result_C_HF, "UHF"; N_up=4, N_down=2)
        save_results("carbon_LDA.txt", bset_C, result_C_LDA, "LSDA"; N_up=4, N_down=2)


        # ============================================================================
        # SUMMARY
        # ============================================================================

        println("\n\n")
        println("="^80)
        println("  SUMMARY OF RESULTS")
        println("="^80)

        println("\nHydrogen Atom (1s¹):")
        @printf "  UHF Total Energy: %.10f Hartree\n" result_H_HF.energy
        println("  (LSDA not applicable for single-electron systems)")

        println("\nCarbon Atom (1s² 2s² 2p², ³P):")
        @printf "  UHF  Total Energy: %.10f Hartree\n" result_C_HF.energy
        @printf "  LSDA Total Energy: %.10f Hartree\n" result_C_LDA.energy

        println("\nDetailed results saved to:")
        println("  - hydrogen_HF.txt")
        println("  - carbon_HF.txt")
        println("  - carbon_LDA.txt")

        println("\n" * "="^80)
        println("  Calculation Complete!")
        println("="^80)

        return 0  # Success
    catch e
        println(stderr, "Error: ", e)
        Base.show_backtrace(stderr, catch_backtrace())
        return 1  # Error
    end
end

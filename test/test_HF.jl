module test_HF

using Test
using DFTatom
using GaussianBasis
using LinearAlgebra

@testset "HF.jl" begin
    @testset "one_electron_energy" begin
        H_bset = BasisSet("sto-3g", "H 0.0 0.0 0.0")
        expected_H_energy = -0.466581 # More precise value for H/sto-3g
        @test DFTatom.one_electron_energy(H_bset) ≈ expected_H_energy atol = 1e-6
    end

    @testset "SCF: Hydrogen" begin
        H_bset = BasisSet("sto-3g", "H 0.0 0.0 0.0")
        results = DFTatom.SCF(H_bset; n_unpaired = 1)
        # For a one-electron system, SCF energy should be the same as one-electron energy
        expected_H_energy = -0.466581
        @test results.energy ≈ expected_H_energy atol = 1e-6
    end

    @testset "SCF: Helium" begin
        He_bset = BasisSet("sto-3g", "He 0.0 0.0 0.0")
        results = DFTatom.SCF(He_bset; n_unpaired = 0)
        # Reference value for He/sto-3g RHF energy is ~-2.85516 Hartree
        expected_He_energy = -2.85516
        @test results.energy ≈ expected_He_energy atol = 1e-5
    end

    @testset "SCF: Carbon" begin
        C_bset = BasisSet("sto-3g", "C 0.0 0.0 0.0")
        results = DFTatom.SCF(C_bset; n_unpaired = 2, tol = 1e-7, maxiter=200)
        # Reference value for C/sto-3g UHF energy (triplet) is ~-37.40666 Hartree
        expected_C_energy = -37.40666
        @test results.energy ≈ expected_C_energy atol = 1e-5
    end
end

end # module test_HF
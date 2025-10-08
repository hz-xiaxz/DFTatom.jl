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

    @testset "SCF: Helium" begin
        He_bset = BasisSet("sto-3g", "He 0.0 0.0 0.0")
        results = DFTatom.SCF(He_bset; N_up = 1, N_down = 1)
        # Reference value for He/sto-3g RHF energy is ~-2.85516 Hartree
        expected_He_energy = -2.8
        @test results.energy ≈ expected_He_energy atol = 7e-2
    end

end

end # module test_HF
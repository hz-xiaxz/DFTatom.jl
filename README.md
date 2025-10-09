# DFTatom

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://hz-xiaxz.github.io/DFTatom.jl/stable)
[![Development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://hz-xiaxz.github.io/DFTatom.jl/dev)
[![Test workflow status](https://github.com/hz-xiaxz/DFTatom.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/hz-xiaxz/DFTatom.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/hz-xiaxz/DFTatom.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/hz-xiaxz/DFTatom.jl)
[![Docs workflow Status](https://github.com/hz-xiaxz/DFTatom.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/hz-xiaxz/DFTatom.jl/actions/workflows/Docs.yml?query=branch%3Amain)
[![BestieTemplate](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JuliaBesties/BestieTemplate.jl/main/docs/src/assets/badge.json)](https://github.com/JuliaBesties/BestieTemplate.jl)

This repository contains a simple implementation of Hartree-Fock (HF) and Density Functional Theory (DFT) for atomic systems, developed as a homework project. The code is written in Julia and demonstrates the self-consistent field (SCF) procedure for calculating atomic energies with various Gaussian basis sets.

# Results

## Carbon Atom (Excited State, ¹D, N_up=3, N_down=3)

Hartree-Fock calculation with a convergence tolerance of 1e-6.

| Basis   | E (a.u.)       | Iterations |
|:--------|:---------------|:-----------|
| sto-3g  | -37.0895866209 | 1          |
| sto-6g  | -37.4635199205 | 1          |
| 6-31g   | -37.7729621437 | 96         |
| 6-31g*  | -37.6817450955 | 96         |
| cc-pvdz | -37.9155567340 | 101        |

# Setup

Before running the scripts, you need to have Julia installed on your system.

1.  **Install Julia**: Follow the instructions on the official Julia website: [julialang.org/install/](https://julialang.org/install/).
2.  **Julia Version**: Please install Julia version **1.11.7**. The project's dependency `GaussianBasis.jl` does not yet support Julia 1.12.

Once Julia is installed, the project's other dependencies will be automatically handled by the `--project` flag when you run the script.

# Reproducing Results

The energy calculations presented in the Results section can be reproduced by running the provided script.

1.  Open your terminal and navigate to the root directory of this project.
2.  Run the following command:

    ```bash
    julia --project scripts/run_hf.jl
    ```

This will execute the Hartree-Fock calculations for the Carbon atom with the different basis sets and print the results table to your console.

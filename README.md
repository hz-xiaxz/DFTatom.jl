# DFTatom

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://hz-xiaxz.github.io/DFTatom.jl/stable)
[![Development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://hz-xiaxz.github.io/DFTatom.jl/dev)
[![Test workflow status](https://github.com/hz-xiaxz/DFTatom.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/hz-xiaxz/DFTatom.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/hz-xiaxz/DFTatom.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/hz-xiaxz/DFTatom.jl)
[![Docs workflow Status](https://github.com/hz-xiaxz/DFTatom.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/hz-xiaxz/DFTatom.jl/actions/workflows/Docs.yml?query=branch%3Amain)
[![BestieTemplate](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JuliaBesties/BestieTemplate.jl/main/docs/src/assets/badge.json)](https://github.com/JuliaBesties/BestieTemplate.jl)

A Julia implementation of **Unrestricted Hartree-Fock (UHF)** and **Density Functional Theory with Local Spin Density Approximation (LSDA)** for atomic systems. This package demonstrates self-consistent field (SCF) calculations for spin-polarized atoms using Gaussian basis sets.

## Features

- **Unrestricted Hartree-Fock (UHF)**: Separate treatment of spin-up and spin-down electrons
- **DFT-LSDA**: Local Spin Density Approximation 
- **Numerical Integration**: Lebedev angular grids (Order 23, 302 points) with logarithmic radial grids (500 points)
- **Symmetry Breaking**: Aufbau-based initial guess with orbital preference for degenerate states
- **Multiple Basis Sets**: Support for STO-3G, STO-6G, 6-31G, 6-31G*, cc-pVDZ
- **Orbital Analysis**: Automatic orbital labeling (1s, 2s, 2p_x, etc.) with wavefunction coefficients

## Results

### Hydrogen Atom (Ground State: 1s¹, ²S)

**UHF Energy**: -0.4982329092 Hartree

*Note: LSDA is not applicable for single-electron systems.*

### Carbon Atom (Ground State: 1s² 2s² 2p², ³P Triplet)

Configuration: N_up=4, N_down=2

#### Hartree-Fock (UHF)

| Basis Set | Energy (Hartree) | Convergence |
|:----------|:-----------------|:------------|
| sto-3g    | -37.1983925466   | 1 iter      |
| sto-6g    | -37.5723640968   | 1 iter      |
| 6-31g*    | -37.6805765603   | 27 iters    |
| cc-pvdz   | -37.6865443947   | 29 iters    |

#### DFT-LDA (LSDA)

| Basis Set | Energy (Hartree) | Convergence |
|:----------|:-----------------|:------------|
| sto-3g    | -36.5852143503   | 1 iter      |
| sto-6g    | -36.9598834496   | 1 iter      |
| 6-31g*    | -37.0939828460   | 539 iters   |
| cc-pvdz   | -37.1027322229   | 297 iters   |

## Installation & Setup

### Prerequisites

Install Julia from the official website: [julialang.org/downloads/](https://julialang.org/downloads/)

### Quick Start

```bash
# Clone the repository
cd /path/to/DFTatom.jl

# Activate the project environment
julia --project=.

# Install dependencies (first time only)
julia --project=. -e "using Pkg; Pkg.instantiate()"

# Run the main assignment script
julia --project=. scripts/assignment.jl
```

## Usage Examples

### Basic Calculations

```julia
using DFTatom
using GaussianBasis

# Hydrogen atom - UHF
bset_H = BasisSet("6-31g_st_", "H 0.0 0.0 0.0")
result_H = SCF(bset_H; N_up=1, N_down=0, maxiter=100, α=0.8)

# Carbon atom - UHF with Aufbau
bset_C = BasisSet("6-31g_st_", "C 0.0 0.0 0.0")
result_C_HF = SCF(bset_C; N_up=4, N_down=2, maxiter=100, α=0.8,
                  use_aufbau=true, favor_high_m=false)

# Carbon atom - LSDA
result_C_LDA = KS_SCF(bset_C; N_up=4, N_down=2, maxiter=1000, α=0.5,
                      use_aufbau=true, favor_high_m=false, tol=1e-5)
```

### Testing Different Basis Sets

```bash
# Compare HF energies across basis sets
julia --project=. scripts/run_hf.jl

# Compare LSDA energies across basis sets
julia --project=. scripts/run_lda.jl
```

## Output Files

Running `scripts/assignment.jl` generates detailed analysis files:
- `hydrogen_HF.txt` - Hydrogen UHF results with orbital energies and wavefunctions
- `carbon_HF.txt` - Carbon UHF results with complete orbital analysis
- `carbon_LDA.txt` - Carbon LSDA results with complete orbital analysis

Each file contains:
- Total energy
- Orbital energies for all occupied orbitals + LUMO
- Orbital labels (1s, 2s, 2p_x, 2p_y, 2p_z, etc.)
- Wavefunction expansion coefficients in the Gaussian basis

## Documentation

For detailed technical information, see [ASSIGNMENT.md](ASSIGNMENT.md), which includes:
- Mathematical formulation of UHF and LSDA
- Numerical integration grid specifications
- Symmetry breaking strategies for degenerate orbitals
- Convergence parameters and troubleshooting
- Physical interpretation of results

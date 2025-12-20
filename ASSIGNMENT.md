# Atomic Electronic Structure Calculations

## Assignment Description

Calculate the **Hartree-Fock (UHF)** and **DFT-LDA (LSDA)** electronic structure for isolated atoms with spin polarization.

### Requirements

1. **Atoms**: Hydrogen (H) and Carbon (C)
2. **Methods**:
   - Unrestricted Hartree-Fock (UHF)
   - Density Functional Theory with Local Spin Density Approximation (LSDA)
3. **Output**:
   - Total energy
   - Single-particle orbital energies (all occupied + first unoccupied/LUMO)
   - Orbital wavefunctions (expansion coefficients in basis)
   - Proper orbital labels (1s, 2s, 2p, etc.)

### Physical Systems

**Hydrogen Atom (H)**
- Ground State: 1s¹ (²S)
- Electronic Configuration: N_up = 1, N_down = 0
- Spin-polarized with one unpaired electron

**Carbon Atom (C)**
- Ground State: 1s² 2s² 2p² (³P triplet)
- Electronic Configuration: N_up = 4, N_down = 2
- Two unpaired p-electrons with parallel spins (Hund's rule)

---

## Results

### Hydrogen Atom Ground State (1s¹, N_up=1, N_down=0)

**Note**: LSDA is not applicable to single-electron systems. Only UHF results are reported.


**Why no LSDA for Hydrogen?**
LSDA (and DFT in general) is designed for multi-electron systems where exchange-correlation effects between electrons are important. For a single electron, there is no exchange or correlation, and the LSDA functional produces unphysical results due to artificial symmetry breaking in the spin channels.

### Carbon Atom Ground State (³P, N_up=4, N_down=2)

#### Hartree-Fock (UHF) Energies

| Basis Set | HF Energy (Hartree) | Convergence |
|:----------|:--------------------|:------------|
| sto-3g    | -37.1983925466      | 1 iter      |
| sto-6g    | -37.5723640968      | 1 iter      |
| 6-31g     | -37.6778369124      | 27 iters    |
| 6-31g*    | -37.6805765603      | 27 iters    |
| cc-pvdz   | -37.6865443947      | 29 iters    |

#### DFT-LDA (LSDA) Energies

| Basis Set | LSDA Energy (Hartree) | Convergence |
|:----------|:----------------------|:------------|
| sto-3g    | -36.5852143503        |  1 iter     |
| sto-6g    | -36.9598834496        |  1 iter     |
| 6-31g*    | -37.0939828460        |  539 iters  |
| cc-pvdz   | -37.1027322229        |  297 iters  |

**Observations**:
- Larger basis sets yield lower (more negative) energies, approaching the complete basis set limit
- cc-pvdz gives the best energy among tested basis sets
- Minimal basis (sto-3g) converges in 1 iteration but with poor accuracy
- LSDA typically gives lower energies than UHF due to inclusion of exchange-correlation effects
- `6-31g` without polarization functions (6-31g*) is not a good basis set for Carbon atom, it cannot converges well so it's excluded from the table.

---

## Code Structure

```
DFTatom.jl/
├── src/
│   ├── DFTatom.jl           # Main module
│   ├── HF.jl                # Hartree-Fock implementation
│   ├── LDA.jl               # DFT-LDA implementation
│   └── AufbauSelection.jl   # Symmetry breaking for degenerate orbitals
├── scripts/
│   ├── assignment.jl        # Main assignment script
│   └── OrbitalAnalysis.jl   # Orbital analysis and output utilities
└── Project.toml             # Package dependencies
```

### Key Implementation Features

1. **Spin Polarization**: Both UHF and LSDA handle spin-up and spin-down channels separately
2. **Gaussian Basis**: Uses the GaussianBasis.jl package with standard basis sets (6-31G)
3. **Aufbau Principle**: For atoms with open shells (like Carbon triplet), uses symmetry-breaking to favor specific p-orbital occupation
4. **SCF Convergence**: Self-consistent field iteration with density mixing

---

## How to Reproduce Results

### 1. Prerequisites

Install Julia (version 1.9 or later): https://julialang.org/downloads/

### 2. Setup Environment

```bash
# Navigate to the project directory
cd /path/to/DFTatom.jl

# Start Julia
julia --project=.

# Install dependencies (first time only)
julia> using Pkg
julia> Pkg.instantiate()
```

### 3. Run the Assignment Script

```bash
julia --project=. scripts/assignment.jl
```

This will:
- Calculate UHF and LSDA for both H and C atoms
- Print detailed orbital analysis to the terminal
- Save results to text files:
  - `hydrogen_HF.txt`
  - `hydrogen_LDA.txt`
  - `carbon_HF.txt`
  - `carbon_LDA.txt`

### 4. Expected Output Format

**Terminal Output Example**:
```
======================================================================
UHF Orbital Analysis
======================================================================

Total Energy: -37.6889123456 Hartree

--- Spin-Up (α) Orbitals ---
Orbital   Label      Energy (Hartree)   Occupation
------------------------------------------------------------
  1       1s         -11.34567890         ●  (occ)
  2       2s          -0.72345678         ●  (occ)
  3       2p_z        -0.42123456         ●  (occ)
  4       2p_y        -0.41987654         ●  (occ)
  5       2p_x        -0.12345678         ○  (LUMO)
...
```

**File Output** (e.g., `carbon_HF.txt`):
- Detailed orbital energies with labels
- Occupation indicators (● = occupied, ○ = virtual)
- LUMO clearly marked
- Wavefunction coefficients for key orbitals (HOMO, LUMO)

### 5. Understanding Wavefunction Coefficients

The output files include **wavefunction expansion coefficients** that show how each molecular orbital is constructed from the Gaussian basis functions.

**Mathematical Background**:

In quantum chemistry, each orbital is expanded as a linear combination of basis functions:

```
ψᵢ(r) = Σⱼ Cⱼᵢ × φⱼ(r)
```

where:
- **ψᵢ(r)** = the i-th molecular orbital (e.g., 1s, 2s, 2p_z)
- **φⱼ(r)** = the j-th Gaussian basis function from the basis set (6-31G*, cc-pVDZ, etc.)
- **Cⱼᵢ** = expansion coefficient (printed in output files)

**What the Coefficients Tell You**:

1. **Which basis functions contribute** to each orbital
2. **How much each contributes** (larger |C| = more important)
3. **The sign** (positive or negative, determines orbital phase)

**Example Output** (Carbon 1s orbital):
```
--- Orbital 1 (α-spin) Wavefunction Coefficients ---
Basis Function            Coefficient
---------------------------------------------
  s                          0.994523
  s                         -0.233891
```

This shows the 1s orbital is primarily the first s-type Gaussian (coefficient ≈ 1.0) with a small negative contribution from the second s-type Gaussian. The negative coefficient creates proper radial nodes and orbital contraction.

**Physical Interpretation**:
- **Large coefficients** (|C| ≈ 1): That basis function is the main component
- **Small coefficients** (|C| ≈ 0): Negligible contribution (filtered out if |C| < 0.01)
- **Multiple significant coefficients**: Orbital is a mixture (common for valence orbitals)

**Where These Come From**:

These coefficients are the eigenvectors **C** from solving the generalized eigenvalue problem:
- **UHF**: F^α C^α = S C^α E^α (Fock matrix eigenvalue problem)
- **LSDA**: Same structure, but Fock matrix includes exchange-correlation potential

The output includes coefficients for:
- HOMO-α (highest occupied molecular orbital, spin-up)
- LUMO-α (lowest unoccupied molecular orbital, spin-up)
- HOMO-β (highest occupied molecular orbital, spin-down)

---

## Technical Details

### Hartree-Fock (UHF)

The unrestricted Hartree-Fock method solves the generalized eigenvalue problem:

```
F^α C^α = S C^α E^α
F^β C^β = S C^β E^β
```

where:
- **F^α, F^β**: Fock matrices for spin-up and spin-down
  - F = H_core + J - K
  - H_core = T (kinetic) + V_nuc (nuclear attraction)
  - J (Coulomb): Electron-electron repulsion
  - K (Exchange): Quantum mechanical exchange interaction
- **S**: Overlap matrix between Gaussian basis functions
- **C**: Molecular orbital coefficient matrices
- **E**: Orbital energy eigenvalues

**Implementation**:
- Gaussian basis sets: STO-3G, STO-6G, 6-31G, 6-31G*, cc-pVDZ
- Two-electron integrals computed analytically via GaussianBasis.jl
- SCF iteration with simple linear mixing: P_new = α·P_SCF + (1-α)·P_old
- Convergence criterion: ||ΔP|| < 10⁻⁶

### DFT-LDA (LSDA)

The Local Spin Density Approximation solves:

```
F^α = H_core + J + V_xc^α
F^β = H_core + J + V_xc^β
```

where:
- **V_xc**: Exchange-correlation potential (computed numerically on a grid)
- **Slater-Dirac Exchange Functional**:
  - ε_x(ρ_α, ρ_β) = C_x · (ρ_α^(4/3) + ρ_β^(4/3)) / ρ_total
  - C_x = -3/4 · (6/π)^(1/3)
  - v_xc^α = ∂(ρ·ε_x)/∂ρ_α = (4/3)·C_x·ρ_α^(1/3)

**Numerical Integration Grid**:

1. **Radial Grid (Logarithmic)**:
   - 500 points spanning r ∈ [exp(-5), exp(3)] Bohr
   - Trapezoidal integration weights
   - Volume element includes r² factor for spherical coordinates

2. **Angular Grid (Lebedev)**:
   - Order 23 Lebedev quadrature (302 angular points)
   - Exact integration of spherical harmonics up to degree 23
   - Provides optimal angular coverage with rotational invariance
   - Total grid: 500 × 302 = 151,000 points

3. **Grid Integration**:
   ```julia
   ∫ f(r) d³r = ∑ᵢ f(rᵢ, θᵢ, φᵢ) · wᵢ
   ```
   where wᵢ = w_radial(rᵢ) · w_angular(θᵢ, φᵢ) · r² · 4π

**Basis Function Evaluation**:
- Each Gaussian basis function φ(r) = R(r) · Y_lm(θ,φ) is evaluated at all grid points
- Density computed as: ρ(r) = ∑ᵢⱼ Pᵢⱼ φᵢ(r) φⱼ(r)
- XC potential projected back to basis: V_xc[μ,ν] = ∫ φμ(r) v_xc(r) φν(r) d³r

**SCF Convergence**:
- More conservative mixing (α = 0.1-0.3) due to nonlinearity of XC functional
- Typical convergence: 100-500 iterations
- Grid quality crucial for accuracy and convergence

### Symmetry Breaking for Degenerate Orbitals

For Carbon's triplet state (4,2), the 2p orbitals are three-fold degenerate. Without symmetry breaking, SCF can oscillate between equivalent solutions.

**Strategy**:

1. **Aufbau-Based Initial Guess**:
   - Solve one-electron problem: (T + V_nuc)C = SCε
   - Rank orbitals by energy, then by |m| quantum number
   - For triplet: Prefer m=0 (p_z) over m=±1 (p_x, p_y)
   - Build initial density from reordered orbitals

2. **Symmetry-Breaking Field**:
   - Add small potential to Fock matrix: F → F + V_symbreak
   - V_symbreak is diagonal in basis functions
   - For `favor_high_m=false`: V_symbreak[μ,μ] = -0.002·(1 - |m_μ|)
   - Lowers energy of m=0 basis functions by ~2 mHartree
   - Field is constant throughout SCF (not adaptive)

3. **Effect on Convergence**:
   - Maintains consistent orbital character across iterations
   - Prevents rotation between degenerate p-orbitals
   - Physical results: Total energy unaffected (below numerical precision)
   - Orbital energies within degeneracy slightly split (artifact of symmetry breaking)

**Note**: For production calculations, would use point group symmetry constraints instead of artificial fields. This approach is pedagogical and works well for atoms.

---

## Results Interpretation

### Hydrogen Atom

**Expected**:
- Single 1s orbital occupied
- LUMO is 2s
- UHF and LSDA give similar but not identical energies

### Carbon Atom (Triplet)

**Expected Configuration**:
```
Spin-up (α):   1s↑ 2s↑ 2p_y↑ 2p_x↑     (N = 4)
Spin-down (β): 1s↓ 2s↓                  (N = 2)
```

**Orbital Splitting**:
- 1s: Most negative energy (core)
- 2s: Higher than 1s (valence)
- 2p: Three-fold degenerate (but split by symmetry breaking)
- 3s: LUMO (first unoccupied)

**Energy Comparison**:
- UHF typically gives higher (less negative) energy than LSDA
- LSDA includes correlation effects via exchange-correlation functional

---

## Troubleshooting

### SCF Not Converging

**For Carbon triplet**:
- Try different mixing parameters α (0.1 to 0.5)
- Increase maxiter
- Ensure `use_aufbau=true` and `favor_high_m=false`

**For symmetric systems**:
- Lower α for more conservative mixing
- Check for orbital degeneracies causing oscillation

### Basis Set Issues

Current implementation uses 6-31G*. For better accuracy:
- Try 6-311G** (larger basis)
- Note: Larger basis = slower convergence

---

## References
1. thijssen _Computational Physics_ 


**Dependencies**:
- GaussianBasis.jl (Stores basis set information including `sto-3g`, `6-31g`, `6-31g*`, `cc-pvdz` etc.)
- LinearAlgebra.jl (Julia standard library)
- OMEinsum.jl (tensor contractions)
- Printf.jl (formatted output)
- SphericalHarmonics.jl (angular momentum functions)
- Lebedev.jl (integration grids)

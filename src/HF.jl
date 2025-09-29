using GaussianBasis
using LinearAlgebra
using OMEinsum

# compute under Hartree unit

"""
    one_electron_energy(bset::BasisSet)

Computes the energy of a single electron in the potential of the nuclei.
This is equivalent to the ground state energy of a one-electron system (like H, He+, Li2+, etc.).

# Arguments
- `bset::BasisSet`: The basis set for the system.

# Returns
- `Float64`: The lowest energy eigenvalue. For a Hydrogen atom with sto-3g, this is ≈ -0.4666 Hartree.
"""
function one_electron_energy(bset::BasisSet)
    T = kinetic(bset)
    nuc = nuclear(bset)
    # For a one-electron system, J and K are zero.
    S = overlap(bset)
    # solving the generalized eigenvalue problem
    evals, evecs = eigen(T + nuc, S)
    return minimum(evals)
end

"""
    init_conf(bset::BasisSet)

Computes an initial guess for the molecular orbital coefficients by solving the one-electron problem (ignoring electron-electron repulsion). The resulting orbitals are sorted by energy.
    
"""
function init_conf(bset::BasisSet)
    # ignore electron repulsion first
    T = kinetic(bset)
    nuc = nuclear(bset)
    H0 = T + nuc
    S = overlap(bset)
    evals, evecs = eigen(H0, S)
    p = sortperm(evals)
    C = evecs[:, p] # For Carbon, px, py, pz orbitals are degenerate
    # return the configuration without interaction
    return C
end

"""
    SCF(bset::BasisSet; n_unpaired::Int, maxiter::Int=100, tol::Float64=1e-6)

Perform a Unrestricted Hartree-Fock (UHF) self-consistent field calculation.
This is suitable for atoms or molecules, including open-shell systems. This implementation assumes the system is neutral (charge=0).

# Arguments
- `bset::BasisSet`: The basis set for the atom.
- `n_unpaired::Int`: The number of unpaired electrons. Must have the same parity as the total number of electrons.
- `maxiter::Int`: Maximum number of SCF iterations.
- `tol::Float64`: Convergence tolerance for the density matrix.

# Returns
- A `NamedTuple` with the following fields:
  - `energy::Float64`: The total SCF energy (electronic + nuclear repulsion).
  - `P_up::Matrix{Float64}`: The alpha-spin density matrix.
  - `P_down::Matrix{Float64}`: The beta-spin density matrix.
  - `C_up::Matrix{Float64}`: The alpha-spin molecular orbital coefficients.
  - `C_down::Matrix{Float64}`: The beta-spin molecular orbital coefficients.
  - `E_up::Vector{Float64}`: The alpha-spin orbital energies.
  - `E_down::Vector{Float64}`: The beta-spin orbital energies.
"""
function SCF(bset::BasisSet; n_unpaired::Int, maxiter = 100, tol = 1e-6)
    # --- Setup ---
    # Assuming a neutral system, charge = 0
    n_electrons = sum(atom.Z for atom in bset.atoms)
    if (n_electrons + n_unpaired) % 2 != 0 || (n_electrons - n_unpaired) < 0
        error("Invalid number of unpaired electrons. n_electrons and n_unpaired must have the same parity, and n_electrons >= n_unpaired.")
    end
    N_up = (n_electrons + n_unpaired) ÷ 2
    N_down = (n_electrons - n_unpaired) ÷ 2

    # Core Hamiltonian (constant)
    T0 = kinetic(bset)
    nuc = nuclear(bset)
    H0 = T0 + nuc
    S = overlap(bset)
    inter = ERI_2e4c(bset)

    # Initial guess orbitals
    C = init_conf(bset)
    C_up = copy(C)
    C_down = copy(C)

    P_up = C_up[:, 1:N_up] * C_up[:, 1:N_up]'
    P_down = C_down[:, 1:N_down] * C_down[:, 1:N_down]'

    # --- SCF iterations ---
    for i = 1:maxiter
        P_total = P_up + P_down

        # Build Fock matrices
        F_up = copy(H0)
        F_down = copy(H0)

        # Coulomb (J)
        @ein! F_up[μ, ν] += P_total[λ, σ] * inter[μ, ν, λ, σ]
        @ein! F_down[μ, ν] += P_total[λ, σ] * inter[μ, ν, λ, σ]

        # Exchange (K) for alpha
        @ein ex_F_up[μ, ν] := P_up[λ, σ] * inter[μ, σ, λ, ν]
        F_up .-= ex_F_up

        # Exchange (K) for beta
        @ein ex_F_down[μ, ν] := P_down[λ, σ] * inter[μ, σ, λ, ν]
        F_down .-= ex_F_down

        # Solve and sort
        evals_up, C_up_new = eigen(F_up, S)
        p_up = sortperm(evals_up)
        evals_up = evals_up[p_up]
        C_up_new = C_up_new[:, p_up]

        evals_down, C_down_new = eigen(F_down, S)
        p_down = sortperm(evals_down)
        evals_down = evals_down[p_down]
        C_down_new = C_down_new[:, p_down]

        # New densities from occupied orbitals
        P_up_new = C_up_new[:, 1:N_up] * C_up_new[:, 1:N_up]'
        P_down_new = C_down_new[:, 1:N_down] * C_down_new[:, 1:N_down]'

        # Convergence
        if norm(P_up_new - P_up) + norm(P_down_new - P_down) < tol
            println("SCF converged in $i iterations.")

            # Energy calculation
            E_elec = 0.5 * (tr((P_up_new + P_down_new) * H0) + tr(P_up_new * F_up) + tr(P_down_new * F_down))
            
            E_nuc = 0.0
            atoms = bset.mol.atoms
            for i_atom in 1:length(atoms)
                for j_atom in (i_atom+1):length(atoms)
                    r_ij = norm(atoms[i_atom].xyz - atoms[j_atom].xyz)
                    if r_ij > 1e-8 # Avoid self-interaction if atoms are at the same position
                        E_nuc += atoms[i_atom].Z * atoms[j_atom].Z / r_ij
                    end
                end
            end
            E_total = E_elec + E_nuc

            return (
                energy = E_total,
                P_up = P_up_new,
                P_down = P_down_new,
                C_up = C_up_new,
                C_down = C_down_new,
                E_up = evals_up,
                E_down = evals_down,
            )
        end

        # Update
        P_up = P_up_new
        P_down = P_down_new
    end

    error("SCF did not converge in $maxiter iterations.")
end
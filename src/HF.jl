using GaussianBasis
using LinearAlgebra
using OMEinsum
using Printf

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
    evals, _ = eigen(T + nuc, S)
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
    SCF(...)

Perform an Unrestricted Hartree-Fock (UHF) self-consistent field calculation.

The UHF equations are solved:
```math
F^α C^α = S C^α E^α \
F^β C^β = S C^β E^β
```
where ``F^α = H_{core} + J - K^α`` and ``F^β = H_{core} + J - K^β`.

# Arguments
- `bset::BasisSet`: The basis set for the atom.
- `N_up::Int`: The number of spin-up electrons.
- `N_down::Int`: The number of spin-down electrons.

# Keyword Arguments
- `maxiter::Int = 100`: Maximum number of SCF iterations.
- `α::Float64 = 0.8`: Mixing parameter for density updates.
- `tol::Float64 = 1e-6`: Convergence tolerance for the density matrix.

# Returns
A tuple containing the total energy, density matrices, and orbital coefficients.
"""
function SCF(
    bset::BasisSet;
    N_up::Int,
    N_down::Int,
    maxiter::Int = 100,
    α::Float64 = 0.8, # mixing parameter
    tol::Float64 = 1e-6,
)

    # Core Hamiltonian (constant)
    T0 = kinetic(bset)
    nuc = nuclear(bset)
    H_core = T0 + nuc
    S = overlap(bset)
    inter = ERI_2e4c(bset)

    # Initial guess orbitals
    C = init_conf(bset)
    C_up = copy(C)
    C_down = copy(C)

    # density matrices
    P_up = C_up[:, 1:N_up] * C_up[:, 1:N_up]'
    P_down = C_down[:, 1:N_down] * C_down[:, 1:N_down]'

    @printf "Starting UHF-SCF calculation with %d basis functions.\n" bset.nbas
    @printf "N_up = %d, N_down = %d\n" N_up N_down

    # --- SCF iterations ---
    for i = 1:maxiter
        P_total = P_up + P_down

        # Build Fock matrices
        F_up = copy(H_core)
        F_down = copy(H_core)

        # Coulomb (J)
        @ein! F_up[μ, ν] += P_total[λ, σ] * inter[μ, ν, λ, σ]
        @ein! F_down[μ, ν] += P_total[λ, σ] * inter[μ, ν, λ, σ]

        # Exchange (K) for alpha
        @ein ex_F_up[μ, ν] := P_up[λ, σ] * inter[μ, λ, σ, ν]
        F_up .-= ex_F_up

        # Exchange (K) for beta
        @ein ex_F_down[μ, ν] := P_down[λ, σ] * inter[μ, λ, σ, ν]
        F_down .-= ex_F_down

        # Solve and sort
        evals_up, C_up_new = eigen(F_up, S)
        evals_down, C_down_new = eigen(F_down, S)

        # New densities from occupied orbitals
        P_up_new = C_up_new[:, 1:N_up] * C_up_new[:, 1:N_up]'
        P_down_new = C_down_new[:, 1:N_down] * C_down_new[:, 1:N_down]'

        # Convergence Check
        if norm(P_up_new - P_up) + norm(P_down_new - P_down) < tol
            @printf "SCF converged in %d iterations.\n" i

            P_total_new = P_up_new + P_down_new
            # Energy calculation
            E_elec =
                0.5 *
                (tr(P_total_new * H_core) + tr(P_up_new * F_up) + tr(P_down_new * F_down))

            return (
                energy = E_elec,
                P_up = P_up_new,
                P_down = P_down_new,
                C_up = C_up_new,
                C_down = C_down_new,
            )
        end

        # DIIS mixing could be implemented here for better convergence
        # For now, simple linear mixing
        P_up = α * P_up_new + (1 - α) * P_up
        P_down = α * P_down_new + (1 - α) * P_down
    end

    error("SCF did not converge in $maxiter iterations.")
end
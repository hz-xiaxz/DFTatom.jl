using GaussianBasis
using LinearAlgebra
using OMEinsum
using Printf
using SphericalHarmonics
using Lebedev

"""
    LSDA_xc(rho_a::Float64, rho_b::Float64)

Computes the exchange energy and potential for a uniform electron gas using the Local Spin Density Approximation (LSDA).
This implementation uses the Slater-Dirac exchange functional.

# Arguments
- `rho_a::Float64`: Spin-up electron density at a point.
- `rho_b::Float64`: Spin-down electron density at a point.

# Returns
- `eps_x`: Exchange energy per particle (ε_x).
- `v_x_a`: Spin-up exchange potential (v_x^α).
- `v_x_b`: Spin-down exchange potential (v_x^β).
"""
function LSDA_xc(rho_a::Float64, rho_b::Float64)
    # A small number to prevent division by zero or NaN for zero density
    TOL = 1e-12
    rho_a = max(rho_a, TOL)
    rho_b = max(rho_b, TOL)

    rho_total = rho_a + rho_b

    # Slater-Dirac exchange energy per particle (epsilon_x)
    Cx = -3.0 / 4.0 * (6.0 / pi)^(1.0 / 3.0)
    eps_x = Cx * (rho_a^(4 / 3) + rho_b^(4 / 3)) / rho_total

    # Functional derivative to get the exchange potential (v_x)
    v_x_a = (4.0 / 3.0) * Cx * rho_a^(1 / 3)
    v_x_b = (4.0 / 3.0) * Cx * rho_b^(1 / 3)

    return eps_x, v_x_a, v_x_b
end

"""
    get_integration_grid(; n_rad::Int=100, n_ang_order::Int=23)

Generates a numerical integration grid for a single atom centered at the origin.
It combines a logarithmic radial grid with a Lebedev angular grid.

# Keyword Arguments
- `n_rad::Int`: Number of radial points.
- `n_ang_order::Int`: The order of the Lebedev angular grid.

# Returns
- `points`: A 3xN array of grid point coordinates.
- `weights`: A vector of N integration weights.
"""
function get_integration_grid(; n_rad::Int = 100, n_ang_order::Int = 23)
    # 1. Radial Grid (logarithmic)
    rad_pts = exp.(range(-5, stop = 3, length = n_rad))
    rad_weights = zeros(n_rad)
    # Use trapezoidal rule for radial weights (dr)
    for i = 2:n_rad-1
        rad_weights[i] = (rad_pts[i+1] - rad_pts[i-1]) / 2
    end
    rad_weights[1] = (rad_pts[2] - rad_pts[1]) / 2
    rad_weights[end] = (rad_pts[end] - rad_pts[end-1]) / 2

    # 2. Angular Grid (Lebedev)
    x_leb, y_leb, z_leb, w_leb = lebedev_by_order(n_ang_order)
    n_ang = length(w_leb)

    # 3. Combine grids
    n_grid = n_rad * n_ang
    grid_points = zeros(3, n_grid)
    grid_weights = zeros(n_grid)

    idx = 1
    for i = 1:n_rad
        r = rad_pts[i]
        # Radial part of the volume element is r^2 * dr
        w_r = r^2 * rad_weights[i]

        for j = 1:n_ang
            # The Lebedev weights are normalized to 1, so the angular part of
            # the volume element is 4π * w_leb
            w_ang = 4π * w_leb[j]

            grid_points[1, idx] = r * x_leb[j]
            grid_points[2, idx] = r * y_leb[j]
            grid_points[3, idx] = r * z_leb[j]

            grid_weights[idx] = w_r * w_ang
            idx += 1
        end
    end

    return grid_points, grid_weights
end

"""
    evaluate_basis(bset::BasisSet, r_vec::AbstractVector{Float64})

Evaluates all basis functions in the `BasisSet` at a given point `r_vec = [x, y, z]`.
This implementation assumes a single atom at the origin.
"""
function evaluate_basis(bset::BasisSet, r_vec::AbstractVector{Float64})
    # For a single atom at the origin, the displacement vector is just the point vector
    x, y, z = r_vec
    r = norm(r_vec)
    θ = (r > 0) ? acos(z / r) : 0.0
    φ = atan(y, x)

    # Pre-compute all spherical harmonics up to the maximum l value in the basis set
    lmax = maximum(shell.l for shell in bset.basis; init = 0)
    Ylm_real = computeYlm(θ, φ; lmax = lmax, SHType = SphericalHarmonics.RealHarmonics())

    bf_values = zeros(bset.nbas)
    bf_idx = 1

    for shell in bset.basis
        l = shell.l
        coeffs = shell.coef
        exps = shell.exp

        # Calculate the radial part of the contracted basis function
        radial_val = 0.0
        for i in eachindex(coeffs)
            radial_val += coeffs[i] * exp(-exps[i] * r * r)
        end

        # Combine with the angular part (solid harmonics)
        for m = -l:l
            # Solid harmonic is r^l * Y_lm
            angular_part = (r^l) * Ylm_real[(l, m)]

            bf_values[bf_idx] = radial_val * angular_part
            bf_idx += 1
        end
    end
    return bf_values
end


"""
    KS_SCF(...)

Perform a Kohn-Sham Self-Consistent Field calculation using LSDA.

# Arguments
- `bset::BasisSet`: The basis set for the atom.
- `N_up::Int`: Number of spin-up electrons.
- `N_down::Int`: Number of spin-down electrons.
- `maxiter::Int`: Maximum number of SCF iterations.
- `α::Float64`: Mixing parameter for density updates.
- `tol::Float64`: Convergence tolerance.

# Returns
A tuple containing the total energy, density matrices, and orbital coefficients.
"""
function KS_SCF(
    bset::BasisSet;
    N_up::Int,
    N_down::Int,
    maxiter::Int = 100,
    α::Float64 = 0.5,
    tol::Float64 = 1e-6,
)
    grid_points, grid_weights = get_integration_grid()
    n_grid = length(grid_weights)
    n_basis = bset.nbas

    T0 = kinetic(bset)
    nuc = nuclear(bset)
    H_core = T0 + nuc
    S = overlap(bset)
    inter = ERI_2e4c(bset)

    C = init_conf(bset)
    C_up = copy(C)
    C_down = copy(C)

    P_up = C_up[:, 1:N_up] * C_up[:, 1:N_up]'
    P_down = C_down[:, 1:N_down] * C_down[:, 1:N_down]'

    phi_on_grid = zeros(n_grid, n_basis)
    for i = 1:n_grid
        phi_on_grid[i, :] = evaluate_basis(bset, grid_points[:, i])
    end

    for iter = 1:maxiter
        P_total = P_up + P_down

        @ein J[μ, ν] := P_total[λ, σ] * inter[μ, ν, λ, σ]

        rho_up_on_grid = sum((phi_on_grid * P_up) .* phi_on_grid, dims = 2)
        rho_down_on_grid = sum((phi_on_grid * P_down) .* phi_on_grid, dims = 2)

        v_xc_up_on_grid = zeros(n_grid)
        v_xc_down_on_grid = zeros(n_grid)
        E_xc = 0.0

        for i = 1:n_grid
            eps_xc, v_xc_up, v_xc_down = LSDA_xc(rho_up_on_grid[i], rho_down_on_grid[i])
            v_xc_up_on_grid[i] = v_xc_up
            v_xc_down_on_grid[i] = v_xc_down
            rho_total_on_grid = rho_up_on_grid[i] + rho_down_on_grid[i]
            E_xc += eps_xc * rho_total_on_grid * grid_weights[i]
        end

        V_xc_up = phi_on_grid' * Diagonal(grid_weights .* v_xc_up_on_grid) * phi_on_grid
        V_xc_down = phi_on_grid' * Diagonal(grid_weights .* v_xc_down_on_grid) * phi_on_grid

        F_up = H_core + J + V_xc_up
        F_down = H_core + J + V_xc_down

        evals_up, C_up_new = eigen(F_up, S)
        evals_down, C_down_new = eigen(F_down, S)

        P_up_new = C_up_new[:, 1:N_up] * C_up_new[:, 1:N_up]'
        P_down_new = C_down_new[:, 1:N_down] * C_down_new[:, 1:N_down]'

        P_total_new = P_up_new + P_down_new
        E_h = 0.5 * tr(P_total_new * J)
        E_one_electron = tr(P_total_new * H_core)
        E_total = E_one_electron + E_h + E_xc
        @show E_total
        if norm(P_up_new - P_up) + norm(P_down_new - P_down) < tol
            @printf "SCF converged in %d iterations.\n" iter
            P_total_new = P_up_new + P_down_new
            E_h = 0.5 * tr(P_total_new * J)
            E_one_electron = tr(P_total_new * H_core)
            E_total = E_one_electron + E_h + E_xc
            return (
                energy = E_total,
                P_up = P_up_new,
                P_down = P_down_new,
                C_up = C_up_new,
                C_down = C_down_new,
                evals_up = evals_up,
                evals_down = evals_down,
            )
        end

        P_up = α * P_up_new + (1 - α) * P_up
        P_down = α * P_down_new + (1 - α) * P_down
    end

    error("SCF did not converge in $maxiter iterations.")
end

"""
    run_lda_h()

Run LSDA calculation for Hydrogen atom.
"""
function run_lda_h()
    bset = BasisSet("sto-3g", "H 0.0 0.0 0.0")
    res = KS_SCF(bset; N_up = 1, N_down = 0)

    homo_energy = res.evals_up[1] # For H, N_up=1, so HOMO is the first spin-up orbital

    @printf "H atom total energy:     %.6f Hartree\n" res.energy
    @printf "H atom 1s orbital energy: %.6f Hartree\n" homo_energy
    return res
end

"""
    run_lda_c()

Run LSDA calculation for Carbon atom.
"""
function run_lda_c()
    bset = BasisSet("sto-3g", "C 0.0 0.0 0.0")
    # Carbon: 1s^2 2s^2 2p^2. N=6. Ground state is triplet (^3P).
    res = KS_SCF(bset; N_up = 4, N_down = 2)
    @printf "C atom energy: %.6f Hartree\n" res.energy
    return res
end

"""
    run_lda_he()

Run LSDA calculation for Helium atom.
"""
function run_lda_he()
    bset = BasisSet("sto-3g", "He 0.0 0.0 0.0")
    # Helium: 1s^2. N=2. Ground state is singlet.
    res = KS_SCF(bset; N_up = 1, N_down = 1)
    @printf "He atom energy: %.6f Hartree\n" res.energy
    return res
end

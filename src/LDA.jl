using GaussianBasis

function LSDA_xc(rho_a::Float64, rho_b::Float64)
    # A small number to prevent division by zero or NaN for zero density
    TOL = 1e-12
    rho_a = max(rho_a, TOL)
    rho_b = max(rho_b, TOL)

    rho_total = rho_a + rho_b

    # Slater-Dirac exchange energy per particle (epsilon_x)
    Cx = -3.0/4.0 * (6.0/pi)^(1.0/3.0)
    eps_x = Cx * (rho_a^(4/3) + rho_b^(4/3)) / rho_total

    # Functional derivative to get the exchange potential (v_x)
    v_x_a = (4.0/3.0) * Cx * rho_a^(1/3)
    v_x_b = (4.0/3.0) * Cx * rho_b^(1/3)

    return eps_x, v_x_a, v_x_b
end

ρ
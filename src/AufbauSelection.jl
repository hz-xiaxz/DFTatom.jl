# Shared functions for Aufbau principle with symmetry breaking
# Used by both HF.jl and LDA.jl

"""
    get_basis_angular_momentum(bset::BasisSet)

Returns the angular momentum quantum numbers (l, m) for each basis function.

# Arguments
- `bset::BasisSet`: The basis set.

# Returns
- `l_values`: Vector of l quantum numbers for each basis function.
- `m_values`: Vector of m quantum numbers for each basis function.
"""
function get_basis_angular_momentum(bset::BasisSet)
    l_values = Int[]
    m_values = Int[]

    for shell in bset.basis
        l = shell.l
        for m in -l:l
            push!(l_values, l)
            push!(m_values, m)
        end
    end

    return l_values, m_values
end


"""
    select_orbitals_aufbau(C::Matrix, evals::Vector, l_values::Vector{Int}, m_values::Vector{Int}, N_occ::Int; favor_high_m::Bool=true)

Select orbitals to occupy using the Aufbau principle with symmetry breaking.
For degenerate p orbitals, this breaks the degeneracy by preferring specific m values.

# Arguments
- `C::Matrix`: Molecular orbital coefficient matrix.
- `evals::Vector`: Orbital energies.
- `l_values::Vector{Int}`: Angular momentum l for each basis function.
- `m_values::Vector{Int}`: Magnetic quantum number m for each basis function.
- `N_occ::Int`: Number of orbitals to occupy.
- `favor_high_m::Bool`: If true, favor high |m| (px/py). If false, favor low |m| (pz).

# Returns
- Vector of indices of orbitals to occupy.

# Details
The function calculates each molecular orbital's |m|-character and uses it to break
degeneracies:
- `favor_high_m=true`: Preferentially occupies orbitals with higher |m| (px/py character),
  leaving m=0 (pz) unoccupied. Useful for singlet states where you want (px)²(py)² configuration.
- `favor_high_m=false`: Preferentially occupies orbitals with lower |m| (pz character),
  leaving px/py unoccupied. Useful for triplet states with specific spatial symmetry.
"""
function select_orbitals_aufbau(C::Matrix, evals::Vector, l_values::Vector{Int}, m_values::Vector{Int}, N_occ::Int; favor_high_m::Bool=true)
    n_orb = length(evals)

    # Calculate the m-character for each molecular orbital
    m_character = zeros(n_orb)
    for i in 1:n_orb
        coeff_sq = C[:, i].^2
        # Weighted average of |m| values (higher means more px/py character)
        m_character[i] = sum(coeff_sq .* abs.(m_values)) / (sum(coeff_sq) + 1e-12)
    end

    # Create a sorting key: (energy, ±|m|_character)
    # favor_high_m=true:  sorts to favor higher |m| = px/py (leaves pz unoccupied)
    # favor_high_m=false: sorts to favor lower |m| = pz (leaves px/py unoccupied)
    if favor_high_m
        sort_keys = [(evals[i], -m_character[i], i) for i in 1:n_orb]
    else
        sort_keys = [(evals[i], m_character[i], i) for i in 1:n_orb]
    end
    sort!(sort_keys)

    # Select the first N_occ orbitals
    selected_indices = [key[3] for key in sort_keys[1:N_occ]]

    return selected_indices
end

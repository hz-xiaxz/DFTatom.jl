# Utility functions for analyzing and printing orbital results

using Printf
using LinearAlgebra
using DFTatom

"""
    assign_orbital_labels(bset::BasisSet, C::Matrix, n_orbitals::Int)

Assign quantum number labels (n, l, m) to molecular orbitals based on their
angular momentum character and energy ordering.

For atoms, orbitals are assigned following the Aufbau filling order:
1s < 2s < 2p < 3s < 3p < 3d < 4s < 4p < 4d < 4f

# Returns
- Vector of tuples (n, l, m_label) for each orbital
"""
function assign_orbital_labels(bset::BasisSet, C::Matrix, n_orbitals::Int)
    labels = []

    # Get basis function quantum numbers
    l_values, m_values = DFTatom.get_basis_angular_momentum(bset)

    # Track how many orbitals of each l type we've assigned (regardless of m)
    l_count = Dict{Int,Int}(0 => 0, 1 => 0, 2 => 0, 3 => 0)

    for i in 1:n_orbitals
        # Determine angular momentum by summing contributions from each l-type
        coeff_sq = C[:, i].^2

        # Calculate total weight for each l value
        l_weights = Dict{Int,Float64}()
        for l in 0:3  # s, p, d, f
            l_weights[l] = sum(coeff_sq[j] for j in 1:length(l_values) if l_values[j] == l; init=0.0)
        end

        # Find dominant l character
        l = argmax(l_weights)

        # Increment count for this l type
        l_count[l] += 1

        # Assign principal quantum number based on how many orbitals of this l we've seen
        # For l=0 (s): 1st is 1s, 2nd is 2s, 3rd is 3s, ...
        # For l=1 (p): 1st-3rd are 2p, 4th-6th are 3p, ...
        # For l=2 (d): 1st-5th are 3d, 6th-10th are 4d, ...
        degeneracy = 2 * l + 1
        shell_num = div(l_count[l] - 1, degeneracy) + 1
        n = l + shell_num

        # For p, d, f orbitals, determine m from the dominant contribution within that l
        m = 0
        if l > 0
            # Find which m component has largest coefficient for this l
            m_weights = Dict{Int,Float64}()
            for m_val in -l:l
                m_weights[m_val] = sum(coeff_sq[j] for j in 1:length(l_values)
                                      if l_values[j] == l && m_values[j] == m_val; init=0.0)
            end
            m = argmax(m_weights)
        end

        # Convert l to letter
        l_letter = l == 0 ? "s" : (l == 1 ? "p" : (l == 2 ? "d" : "f"))

        # Add m subscript for p, d, f orbitals
        if l == 0
            m_label = ""
        elseif l == 1
            m_label = m == -1 ? "y" : (m == 0 ? "z" : "x")
        else
            m_label = string(m)
        end

        push!(labels, (n, l_letter, m_label))
    end

    return labels
end


"""
    print_orbital_analysis(bset::BasisSet, result, method_name::String)

Print detailed orbital analysis including energies, labels, and occupations.

# Arguments
- `bset::BasisSet`: The basis set
- `result`: SCF or KS_SCF result tuple
- `method_name::String`: Name of the method (e.g., "UHF", "LSDA")
"""
function print_orbital_analysis(bset::BasisSet, result, method_name::String; N_up::Int, N_down::Int)
    println("\n" * "="^70)
    println("$method_name Orbital Analysis")
    println("="^70)

    # Get orbital labels
    n_orbitals = length(result.evals_up)
    labels_up = assign_orbital_labels(bset, result.C_up, n_orbitals)
    labels_down = assign_orbital_labels(bset, result.C_down, n_orbitals)

    # Print total energy
    @printf "\nTotal Energy: %.10f Hartree\n" result.energy

    # Print spin-up orbitals
    println("\n--- Spin-Up (α) Orbitals ---")
    println("Orbital   Label      Energy (Hartree)   Occupation")
    println("-" * "^"^60)

    for i in 1:min(N_up + 3, n_orbitals)  # Show occupied + a few unoccupied
        n, l, m = labels_up[i]
        label = m == "" ? "$(n)$(l)" : "$(n)$(l)_$(m)"
        occ = i <= N_up ? "●" : "○"
        status = i <= N_up ? "occ" : (i == N_up + 1 ? "LUMO" : "virt")

        @printf "%3d       %-8s   %16.8f      %s  (%s)\n" i label result.evals_up[i] occ status
    end

    # Print spin-down orbitals
    println("\n--- Spin-Down (β) Orbitals ---")
    println("Orbital   Label      Energy (Hartree)   Occupation")
    println("-" * "^"^60)

    for i in 1:min(N_down + 3, n_orbitals)  # Show occupied + a few unoccupied
        n, l, m = labels_down[i]
        label = m == "" ? "$(n)$(l)" : "$(n)$(l)_$(m)"
        occ = i <= N_down ? "●" : "○"
        status = i <= N_down ? "occ" : (i == N_down + 1 ? "LUMO" : "virt")

        @printf "%3d       %-8s   %16.8f      %s  (%s)\n" i label result.evals_down[i] occ status
    end

    println("\n" * "="^70)
end


"""
    print_wavefunction_coefficients(bset::BasisSet, C::Matrix, orbital_idx::Int, spin::String)

Print the wavefunction expansion coefficients for a specific orbital.

# Arguments
- `bset::BasisSet`: The basis set
- `C::Matrix`: Orbital coefficient matrix
- `orbital_idx::Int`: Which orbital to print
- `spin::String`: "α" or "β"
"""
function print_wavefunction_coefficients(bset::BasisSet, C::Matrix, orbital_idx::Int, spin::String)
    println("\n--- Orbital $orbital_idx ($spin-spin) Wavefunction Coefficients ---")
    println("Basis Function            Coefficient")
    println("-" * "^"^45)

    l_values, m_values = get_basis_angular_momentum(bset)

    bf_idx = 1
    for (shell_idx, shell) in enumerate(bset.basis)
        l = shell.l
        l_letter = l == 0 ? "s" : (l == 1 ? "p" : (l == 2 ? "d" : "f"))

        for m in -l:l
            m_label = l == 0 ? "" : (l == 1 ? (m == -1 ? "y" : (m == 0 ? "z" : "x")) : string(m))
            label = m_label == "" ? "$(l_letter)" : "$(l_letter)_$(m_label)"

            coeff = C[bf_idx, orbital_idx]
            if abs(coeff) > 0.01  # Only print significant coefficients
                @printf "  %-20s   %10.6f\n" label coeff
            end

            bf_idx += 1
        end
    end
end


"""
    save_results(filename::String, bset::BasisSet, result, method_name::String; N_up::Int, N_down::Int)

Save orbital analysis to a file.
"""
function save_results(filename::String, bset::BasisSet, result, method_name::String; N_up::Int, N_down::Int)
    open(filename, "w") do io
        # Redirect stdout to file
        old_stdout = stdout
        redirect_stdout(io)

        print_orbital_analysis(bset, result, method_name; N_up=N_up, N_down=N_down)

        # Print a few key wavefunctions
        println("\n\nKey Orbital Wavefunctions:")
        println("="^70)

        # Print HOMO-α
        print_wavefunction_coefficients(bset, result.C_up, N_up, "α")

        # Print LUMO-α if exists
        if N_up < size(result.C_up, 2)
            print_wavefunction_coefficients(bset, result.C_up, N_up + 1, "α")
        end

        # Print HOMO-β
        if N_down > 0
            print_wavefunction_coefficients(bset, result.C_down, N_down, "β")
        end

        redirect_stdout(old_stdout)
    end

    println("Results saved to $filename")
end

using GaussianBasis

# compute under Hartree unit

H_bset = BasisSet("sto-3g", "H 0.0 0.0 0.0")
C_bset = BasisSet("sto-3g", "C 0.0 0.0 0.0")

function hf_energy(bset::BasisSet)
    # T = kinetic energy integrals
    # T = <phi|-1/2 ∇^2|phi>
    T = kinetic(bset, bset)[1, 1, 1, 1]
    # the hartree and Fock terms are incorporated in the two-electron integrals
    V = ERI_2e4c(bset) # for H it produces a 1x1x1x1 array
    nuc = nuclear(bset)[1, 1, 1, 1]
    J = 0
    K = 0
    if size(V) == (1, 1, 1, 1)
        # for H, the two-electron integral is just the Coulomb repulsion
        J = 0
        K = 0
    else
        for (indices, val) in enumerate(V)
            if indices[1] == indices[2] && indices[3] == indices[4]
                J += val
            elseif indices[1] == indices[4] &&
                   indices[2] == indices[3] &&
                   indices[1] != indices[2]
                K += val
            end
        end
    end

    @show T, J, K, nuc
    energy = T + J - K + nuc
    return energy
end
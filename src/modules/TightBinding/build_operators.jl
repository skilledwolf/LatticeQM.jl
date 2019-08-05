
function get_projector(lat::Lattice, name::String)
    if name == "spin"
        return get_spin_projector(lat)
    else
        error("Requested projector not defined.")
    end

end

function get_spin_projector(lat::Lattice)
    N = atom_count(lat)

    f(k, ψ, ϵ) = real.(ψ' * kron(σZ, Diagonal(ones(N))) * ψ)

    f
end



function M_n(lat::Lattice, n::Vector{Float64})
    N = atom_count(lat)

    kron( sum(n[i] .* σs[i] for i=1:3) , Diagonal(ones(N)))
end

Mx(lat::Lattice) = M_n(lat, [1.0, 0.0, 0.0])
My(lat::Lattice) = M_n(lat, [0.0, 1.0, 0.0])
Mz(lat::Lattice) = M_n(lat, [0.0, 0.0, 1.0])

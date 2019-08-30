## Save data
function exportdata(filename::String, lat::Lattice)
    h5open(filename, isfile(filename) ? "r+" : "w") do file
        g = g_create(file, "Lattice")
        g["A"] = lat.A
        g["atoms"] = lat.atoms
        g["atoms_aux"] = lat.atoms_aux

        g2 = g_create(g, "extradimensions") # create a group
        g2["names"] = collect(keys(lat.extradimensions))              # create a scalar dataset inside the group
        g2["indices"] = collect(values(lat.extradimensions))
    end
end

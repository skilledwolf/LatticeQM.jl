
using ..Structure: Lattice, repeat_atoms, get_filtered_positions

export unitcellplot_layers, unitcellplot_layers!


function unitcellplot_layers(lat0::Lattice; repeat=[0:1,0:1])

    p = plot()
    unitcellplot_layers!(p, lat0, repeat=repeat)

    return p
end

function unitcellplot_layers!(p, lat0::Lattice; repeat=[0:1,0:1])
    lat = deepcopy(lat0)
    repeat_atoms!(lat, repeat)

    # Plot the layers
    atoms1 = get_filtered_positions(lat, "layer", x->x==0.0)
    atoms2 = get_filtered_positions(lat, "layer", x->x==1.0)

    scatter!(p, atoms1[1,:], atoms1[2,:], color=:darkblue, aspect_ratio=:equal)
    scatter!(p, atoms2[1,:], atoms2[2,:], color=:cornflowerblue, aspect_ratio=:equal)

    return p
end

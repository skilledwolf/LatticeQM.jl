
h5myopen(filename::String) = h5open(filename, isfile(filename) ? "r+" : "w")

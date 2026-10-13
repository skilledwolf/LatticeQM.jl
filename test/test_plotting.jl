using Test
using LatticeQM
import LatticeQM.Plotting: RecipesBase

# Exercise recipes directly; no plotting backend is needed for data invariants.
RecipesBase.is_key_supported(::Symbol) = true

@testset "Plotting preserves band data and observable alignment" begin
    path = kpath(Geometries.honeycomb(); num_points=6)
    bands = reshape(collect(1.0:18.0), 3, 6)
    obs = cat(reverse(bands; dims=1), 2 .* bands; dims=3)
    data = Spectrum.BandData(copy(bands), copy(obs), path)
    for n in 1:2
        recipe = only(RecipesBase.apply_recipe(Dict{Symbol,Any}(), data, n))
        plotted = transpose(only(recipe.args))
        colors = transpose(recipe.plotattributes[:marker_z])
        for k in axes(bands, 2)
            @test sort(collect(zip(plotted[:, k], colors[:, k]))) ==
                  sort(collect(zip(bands[:, k], obs[:, k, n])))
        end
        @test data.bands == bands
        @test data.obs == obs
    end
    constant = Spectrum.BandData(bands, ones(3, 6, 1), path)
    recipe = only(RecipesBase.apply_recipe(Dict{Symbol,Any}(:sharpen => 2.0), constant))
    @test all(isfinite, recipe.plotattributes[:marker_z])
end

@testset "Plotting repeats one color per orbital" begin
    lat = Geometries.honeycomb()
    colors = [0.2, 0.8]
    for (cells, copies) in ((0:0, 1), (0:1, 4))
        recipe = only(RecipesBase.apply_recipe(Dict{Symbol,Any}(:supercell => cells), lat, colors))
        @test recipe.plotattributes[:marker_z] == repeat(colors, copies)
        @test length(recipe.args[1]) == length(recipe.plotattributes[:marker_z])
    end
    recipe = only(RecipesBase.apply_recipe(Dict{Symbol,Any}(:sort => "sublattice", :filter => ("sublattice", iszero)), lat))
    @test length(recipe.args[1]) == 1
    @test Lattices.countorbitals(lat) == 2
end

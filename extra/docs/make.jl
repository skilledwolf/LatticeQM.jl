push!(LOAD_PATH,"../../src/")
using Documenter, LatticeQM

ENV["GKSwstype"] = "100" # allows plots without window server

makedocs(
	sitename="LatticeQM.jl",
	format = Documenter.HTML(
        	prettyurls = get(ENV, "CI", nothing) == "true"
	),
	pages = [
        "Home" => "index.md",
        "Getting started" => "getting_started.md"
    ]
)

push!(LOAD_PATH,"../../src/")
using Documenter
using Documenter.Remotes
using LatticeQM

DocMeta.setdocmeta!(LatticeQM, :DocTestSetup, :(using LatticeQM, Plots, Random, LinearAlgebra))

ENV["GKSwstype"] = "100" # allows plots without window server

makedocs(
    sitename = "LatticeQM.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        edit_link = "main",
        inventory_version = "dev"
    ),
    modules = [LatticeQM],
    checkdocs = :exports,
    warnonly = [:missing_docs],
    repo = Remotes.GitLab("skilledwolf", "LatticeQM.jl"),
    pages = Any[
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Tutorials" => Any[
            "Overview" => "tutorials/index.md",
            "Tutorial 1: Structure & Geometry" => "tutorials/structure.md",
            "Tutorial 2: Band Structures" => "tutorials/bands.md",
            "Tutorial 3: Haldane Model & Topology" => "tutorials/hall_effect.md",
            "Tutorial 4: Twisted Bilayer Graphene" => "tutorials/twisted_bilayers.md",
            "Tutorial 5: Hofstadter Butterfly" => "tutorials/hofstadter.md",
            "Tutorial 6: Mean-field Self-Consistency" => "tutorials/meanfield.md",
            "Tutorial 7: Floquet Dynamics" => "tutorials/floquet.md",
            "Tutorial 8: Superconductivity (BdG)" => "tutorials/superconductivity.md"
        ],
        "Concept Guides" => Any[
            "Geometry & Lattices" => "guides/geometry.md",
            "Operators & Hamiltonians" => "guides/operators.md",
            "Observables & Spectral Analysis" => "guides/observables.md",
            "Mean-field & Superconductivity" => "guides/meanfield.md",
            "Floquet Workflows" => "guides/floquet.md"
        ],
        "Advanced Workflows" => "advanced/index.md",
        "Contributor Guide" => "contributing.md",
        "API Reference" => "api.md",
        "Publications" => "publications.md"
    ]
)

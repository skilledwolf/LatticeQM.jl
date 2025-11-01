# API Reference

This reference is generated from in‑source docstrings. Use Julia help mode (`?name`)
in the REPL for the most up‑to‑date signatures and keyword defaults. Where
available, cross‑links to tutorials and guides provide usage context.

```@meta
CurrentModule = LatticeQM
```

The sections below list public and internal names grouped by submodule. Use Julia’s help
mode (`?`) for inline documentation and follow cross-links back to tutorials and
guides for practical context.

## Structure
```@autodocs
Modules = [LatticeQM.Structure, LatticeQM.Structure.Lattices, LatticeQM.Structure.Geometries]
Order = [:type, :function, :macro]
Private = true
```

```@docs
LatticeQM
LatticeQM.Structure
LatticeQM.Structure.Lattices
LatticeQM.Structure.Geometries
```

## TightBinding
```@autodocs
Modules = [LatticeQM.TightBinding]
Order = [:type, :function]
Private = true
```

## Operators
```@autodocs
Modules = [LatticeQM.Operators]
Order = [:type, :function]
Private = true
```

## Spectrum
```@autodocs
Modules = [LatticeQM.Spectrum]
Order = [:type, :function]
Private = true
```

## LinearResponse
```@autodocs
Modules = [LatticeQM.LinearResponse]
Order = [:type, :function]
Private = true
```

## Meanfield
```@autodocs
Modules = [LatticeQM.Meanfield]
Order = [:type, :function]
Private = true
```

## Superconductivity
```@autodocs
Modules = [LatticeQM.Superconductivity]
Order = [:type, :function]
Private = true
```

## Floquet
```@autodocs
Modules = [LatticeQM.Floquet]
Order = [:type, :function]
Private = true
```

## Plotting
```@autodocs
Modules = [LatticeQM.Plotting]
Order = [:type, :function]
Private = true
```

## Utils
```@autodocs
Modules = [LatticeQM.Utils]
Order = [:type, :function]
Private = true
```

```@docs
LatticeQM.Utils.@scalar2vector
```

## Common Entrypoints

```@docs
LatticeQM.getbands
LatticeQM.Spectrum.getdos
LatticeQM.LinearResponse.opticalconductivity
LatticeQM.Meanfield.solvehartreefock
LatticeQM.Meanfield.solveselfconsistent
LatticeQM.Operators.setfilling!
```

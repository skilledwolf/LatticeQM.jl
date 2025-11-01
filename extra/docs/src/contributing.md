# Contributor Guide

Welcome to the development side of LatticeQM. This guide outlines how to
contribute code, documentation, and examples while keeping the project
maintainable.

## Development setup
1. Fork/clone the repository and develop the package locally:
   ```julia
   julia
   using Pkg
   Pkg.develop(path="/path/to/LatticeQM")
   ```
2. Install optional tooling: `Revise`, `Documenter`, `Plots`, `IJulia`.
3. Run the test suite (`julia --project -e 'using Pkg; Pkg.test()'`) before
   opening a merge request.

## Coding standards
- Prefer explicit imports (`import .Module: function`) and keep docstrings up to
  date with argument/keyword descriptions.
- Use ASCII unless physics notation requires otherwise.
- Add concise comments ahead of non-trivial code blocks explaining the intent.

## Documentation workflow
- Update or create doc pages in `extra/docs/src/`.
- When adding notebooks or examples, include a short README summarising usage
  and runtime expectations.
- Run `julia --project=extra/docs/ make.jl` to preview the docs locally.

## Branching and commits
- Work on topic branches (`docs/tutorial-structure`, `feat/meanfield-update`).
- Keep commits focused; describe both intent and scope.
- Reference related issues in commit messages when relevant.

## Review checklist
- Tests pass and new functionality has coverage (unit tests or runnable
  examples).
- Documentation reflects API changes and links to newly created artefacts.
- Examples run end-to-end or include notes about outstanding validation steps.

## Release coordination
- Before tagging a release, ensure the docs build cleanly and all relevant
  hyperlinks resolve.

## Getting help
- Open a discussion or issue with context (Julia version, environment, stack
  trace).
- Tag module maintainers when code touches their area.

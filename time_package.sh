#! /usr/bin/env bash

julia -e 'import Pkg; Pkg.precompile()'
julia -e '@time import LatticeQM'

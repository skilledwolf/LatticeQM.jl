To run on TACC's slurm cluster (at time of writing lonestar6):

1) Install julia-1.10 or newer using `juliaup`.
2) At time of writing, make sure pcre2 module is loaded `module load pcre2`.
3) `module load gcc/13` to avoid issues with linker ENV flags of intel 
4) Note that on lonestar6, even to install and precompile packages one needs to run on the node, otherwise there are not enough ressources
5) Generally, it is a good idea to use `JULIA_DEPOT_PATH` and `JULIAUP_DEPOT_PATH` to manage where the julia binaries and files are stored.
6) Don't forget to set `OMP_NUM_THREADS`

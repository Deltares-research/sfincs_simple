#!/bin/sh
MANPATH=$MANPATH:/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/compilers/man; export MANPATH
PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/compilers/bin:$PATH; export PATH

module load nvidia

#nvfortran -acc=gpu -Minfo=accel -fast -O3 -gpu=ccall sfincs_unstructured/sfincs_unstructured.f90 -o unstructured.exe
nvfortran -acc=gpu -Minfo=accel -fast -O3 -gpu=ccall sfincs_structured/sfincs_structured.f90 -o structured.exe

# input arguments: nmax, mmax, nt, copy, momentum, continuity

#./unstructured.exe 4096 4096 1000 T T T
./structured.exe 4096 4096 1000 T T T

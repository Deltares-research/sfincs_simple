#!/bin/sh
MANPATH=$MANPATH:/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/compilers/man; export MANPATH
PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/compilers/bin:$PATH; export PATH

#nvfortran -acc=gpu -Minfo=accel -fast -O3 -gpu=ccall sfincs_unstructured/sfincs_unstructured.f90 -o unstructured.exe
#nvfortran -acc=gpu -Minfo=accel -fast -O3 -gpu=ccall sfincs_structured/sfincs_structured.f90 -o structured.exe

nvfortran -Minfo=accel -fast -O3 sfincs_structured/sfincs_structured.f90 -o structured_cpu.exe

# input arguments: nmax, mmax, nt, copy, momentum, continuity

#./unstructured.exe 1000 1000 1000 T T T
#./structured.exe 1000 1000 1000 T T T
./structured_cpu.exe 1000 1000 1000 T T T

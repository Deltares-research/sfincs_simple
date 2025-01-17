#!/bin/sh

#SBATCH  -p gpu-shared
#SBATCH  --gres=gpu:1
#SBATCH  --job-name=SFINCS_test

module load nvidia

nvfortran -mp=gpu -Minfo=mp -fast sfincs_structured.f90 -o sfincs_structured_OMP -O3
nvfortran -acc=gpu -Minfo=accel -fast sfincs_structured.f90 -o sfincs_structured_ACC -O3
nvfortran -fast sfincs_structured.cuf -o sfincs_structured_cuda -O3

./sfincs_structured_OMP 64 64 1000 T T T     > OMP_0064.txt
./sfincs_structured_OMP 128 128 1000 T T T   > OMP_0128.txt
./sfincs_structured_OMP 256 256 1000 T T T   > OMP_0256.txt
./sfincs_structured_OMP 512 512 1000 T T T   > OMP_0512.txt
./sfincs_structured_OMP 1024 8192 1000 T T T > OMP_1024.txt
./sfincs_structured_OMP 2048 8192 1000 T T T > OMP_2048.txt
./sfincs_structured_OMP 4096 8192 1000 T T T > OMP_4096.txt
./sfincs_structured_OMP 8192 8192 1000 T T T > OMP_8192.txt

./sfincs_structured_ACC 64 64 1000 T T T     > ACC_0064.txt
./sfincs_structured_ACC 128 128 1000 T T T   > ACC_0128.txt
./sfincs_structured_ACC 256 256 1000 T T T   > ACC_0256.txt
./sfincs_structured_ACC 512 512 1000 T T T   > ACC_0512.txt
./sfincs_structured_ACC 1024 8192 1000 T T T > ACC_1024.txt
./sfincs_structured_ACC 2048 8192 1000 T T T > ACC_2048.txt
./sfincs_structured_ACC 4096 8192 1000 T T T > ACC_4096.txt
./sfincs_structured_ACC 8192 8192 1000 T T T > ACC_8192.txt

./sfincs_structured_cuda 64 64 1000 T T T     > cuda_0064.txt
./sfincs_structured_cuda 128 128 1000 T T T   > cuda_0128.txt
./sfincs_structured_cuda 256 256 1000 T T T   > cuda_0256.txt
./sfincs_structured_cuda 512 512 1000 T T T   > cuda_0512.txt
./sfincs_structured_cuda 1024 8192 1000 T T T > cuda_1024.txt
./sfincs_structured_cuda 2048 8192 1000 T T T > cuda_2048.txt
./sfincs_structured_cuda 4096 8192 1000 T T T > cuda_4096.txt
./sfincs_structured_cuda 8192 8192 1000 T T T > cuda_8192.txt


#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=12:00:00
#SBATCH --partition=gpu_devel
#SBATCH --mem-per-cpu=10G
#SBATCH --gpus=1 
#SBATCH -J  vampire
#SBATCH -o  output
#SBATCH -e  error

# module load apps/vampire/5.0.0; mpirun -np 1 vampire-parallel
module load apps/vampire/5.0.0; vampire-serial
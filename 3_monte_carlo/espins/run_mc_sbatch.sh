#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=06:00:00
#SBATCH --partition=gpu_devel
#SBATCH --mem-per-cpu=10G
#SBATCH --gpus=1 
#SBATCH -J  mc.x
#SBATCH -o  output
#SBATCH -e  error

module load compilers/intel; mpirun -np 1 /trinity/home/d.poletaev/soft/ESpinS/mc.x Fe2O3_afm1
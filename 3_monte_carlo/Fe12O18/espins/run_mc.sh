#!/bin/bash
#SBATCH -p lenovo
#SBATCH -t 1-00:00:00 
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -J  espins
#SBATCH -e  error

module load intel2021/mkl/latest; export UCX_TLS=ud,sm,self; /home/dpoletaev/soft/ESpinS/mc.x Fe12O18_afm1
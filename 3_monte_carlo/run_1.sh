#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=01:00:00
#SBATCH --partition=gpu_devel
#SBATCH --mem-per-cpu=10G
#SBATCH --gpus=1 
#SBATCH -J  1_coupling_constants.py
#SBATCH -o  output
#SBATCH -e  error

source /trinity/shared/opt/anaconda3/bin/activate automag
python 1_coupling_constants.py
#!/usr/bin/bash

#SBATCH --job-name=vento_cfd
#SBATCH --partition=tupi
#SBATCH --nodes=1
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --gres=gpu:1

set -e

mkdir -p sw hw
nvfortran --version > ./sw/nvfortran.txt
env > ./sw/env.txt

nvidia-smi > ./hw/gpu.txt
lscpu > ./hw/cpu.txt
free > ./hw/memory.txt

make cfd 
./cfd

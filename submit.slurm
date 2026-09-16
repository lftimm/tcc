#!/usr/bin/bash

#SBATCH --job-name=vento_cfd
#SBATCH --partition=tupi
#SBATCH --nodes=1
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --gres=gpu:1

make cfd
./cfd
NUM_ARQUIVOS=grep NTR INTERFACE.dat | awk '{print $1}'
make vtk && ./build_vtks.sh tcc_lt 1 $NUM_ARQUIVOS

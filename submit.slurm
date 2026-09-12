#!/usr/bin/bash

#SBATCH --job-name=vento_cfd
#SBATCH --partition=tupi
#SBATCH --nodes=1
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=24:00:00

make cuda
./cfd

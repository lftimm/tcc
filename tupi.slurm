#!/usr/bin/bash

#SBATCH --job-name=cfd_analise_desempenho
#SBATCH --partition=tupi
#SBATCH --nodes=1
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=3:00:00
#SBATCH --gres=gpu:1

export PATH=$NVHPC_ROOT/compilers/bin:$PATH

./run.sh on_ssh

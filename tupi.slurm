#!/usr/bin/bash

#SBATCH --job-name=vento_cfd
#SBATCH --partition=tupi
#SBATCH --nodes=1
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=15:00:00
#SBATCH --gres=gpu:1

set -e

export PATH=$NVHPC_ROOT/compilers/bin:$PATH

mkdir -p sw hw

nvc --version > ./sw/nvc.txt
nvfortran --version > ./sw/nvfortran.txt
env > ./sw/env.txt
echo "NVHPC_ROOT=$NVHPC_ROOT" > ./sw/nvhpc_root.txt

nvidia-smi > ./hw/gpu.txt
nvidia-smi -q > ./hw/gpu_full_query.txt
nvidia-smi --query-gpu=name,driver_version,compute_cap,memory.total --format=csv > ./hw/gpu_identity.csv
lscpu > ./hw/cpu.txt
free > ./hw/memory.txt

QUERY="timestamp,utilization.gpu,utilization.memory,memory.used,memory.total"
nvidia-smi --query-gpu=$QUERY --format=csv -l 1 > ./sw/nvidia-smi.csv &
NVSMI_PID=$!

nvidia-smi --query-compute-apps=timestamp,pid,process_name,used_memory --format=csv -l 1 > ./sw/nvidia-smi-procs.csv &
NVSMI_PROCS_PID=$!

cleanup() {
  kill "$NVSMI_PID" "$NVSMI_PROCS_PID" 2>/dev/null || true
}
trap cleanup EXIT

make cfd > ./sw/compilation.txt

# Record which GPU architecture(s) actually got compiled into the binary - if
# the running GPU's compute capability isn't in this list, the driver is
# JIT-recompiling from PTX at every kernel launch instead of running native
# code, which depends on the NVIDIA HPC SDK version used to build.
{ cuobjdump --list-elf ./cfd 2>&1 || echo "cuobjdump not available"; } > ./sw/cfd_embedded_sass.txt
{ cuobjdump --list-ptx ./cfd 2>&1 || echo "cuobjdump not available"; } > ./sw/cfd_embedded_ptx.txt

set +e
./cfd
CFD_EXIT=$?
set -e

nvidia-smi -q > ./hw/gpu_full_query_after.txt
{ dmesg -T 2>&1 | tail -n 200 || echo "dmesg not permitted in this environment"; } > ./sw/dmesg_tail.txt

echo "$CFD_EXIT" > ./sw/cfd_exit_code.txt
if [ "$CFD_EXIT" -eq 0 ]; then
  echo "cfd finished normally."
else
  echo "cfd exited with code $CFD_EXIT - check ./sw/ and ./hw/ for diagnostics."
fi

exit "$CFD_EXIT"

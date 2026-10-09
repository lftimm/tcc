#!/usr/bin/bash

set -e

if [ -n "$NVHPC_ROOT" ]; then
  export PATH="$NVHPC_ROOT/compilers/bin:$PATH"
fi

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


set +e
./cfd > ./sw/out.txt
CFD_EXIT=$?
set -e

nvidia-smi -q > ./hw/gpu_full_query_after.txt

echo "$CFD_EXIT" > ./sw/cfd_exit_code.txt
if [ "$CFD_EXIT" -eq 0 ]; then
  echo "conclusão: cfd terminou de rodar."
  tail ./sw/out.txt
else
  echo "erro: $CFD_EXIT - ./sw/ e ./hw/ tem informações de diagnóstico."
fi

exit "$CFD_EXIT"

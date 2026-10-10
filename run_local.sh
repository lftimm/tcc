#!/usr/bin/bash

set -e

if [ -n "$NVHPC_ROOT" ]; then
  export PATH="$NVHPC_ROOT/compilers/bin:$PATH"
fi

echo "[info] Criando pastas com informações de software e hardware"
rm -rf sw hw
mkdir -p sw hw

nvc --version > ./sw/nvc.txt
nvfortran --version > ./sw/nvfortran.txt
env > ./sw/env.txt
echo "NVHPC_ROOT=$NVHPC_ROOT" > ./sw/nvhpc_root.txt

echo "[info] Coletando informações sobre o ambiente cuda"
nvidia-smi > ./hw/gpu.txt
nvidia-smi -q > ./hw/gpu_full_query.txt
nvidia-smi --query-gpu=name,driver_version,compute_cap,memory.total --format=csv > ./hw/gpu_identity.csv

echo "[info] Coletando informações de memória e cpu"
lscpu > ./hw/cpu.txt
free > ./hw/memory.txt

echo "[info] Utilizando nvidia-smi para pegar informações de uso de gpu"
QUERY="timestamp,utilization.gpu,utilization.memory,memory.used,memory.total"
nvidia-smi --query-gpu=$QUERY --format=csv -l 1 > ./sw/nvidia-smi.csv &
NVSMI_PID=$!

cleanup() {
  kill "$NVSMI_PID" 2>/dev/null || true
}
trap cleanup EXIT

echo "[run] Compilando código cuda"
make cfd > ./sw/compilation.txt

set +e
echo "[run] Rodando cfd, o output ficará armazenado em ./sw/out.txt"
./cfd | tee -a ./sw/out.txt
CFD_EXIT=$?
set -e

echo "[info] Coletando informações após o uso do código"
nvidia-smi -q > ./hw/gpu_full_query_after.txt

echo "$CFD_EXIT" > ./sw/cfd_exit_code.txt
if [ "$CFD_EXIT" -eq 0 ]; then
  echo "[end] Cfd terminou de rodar com êxito, verificar pastas sw e hw para mais informações da run."
  tail ./sw/out.txt
else
  echo "[end] Erro de cfd, $CFD_EXIT - ./sw/ e ./hw/ tem informações de diagnóstico."
fi

exit "$CFD_EXIT"


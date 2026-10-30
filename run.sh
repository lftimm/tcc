#!/usr/bin/bash

set -e

ON_SSH=0
PROFILE=0
for arg in "$@"; do
  case "$arg" in
    on_ssh) ON_SSH=1 ;;
    profile) PROFILE=1 ;;
  esac
done

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

echo "[info] Registrando bibliotecas dinamicamente ligadas ao executavel"
ldd ./cfd > ./sw/cfd_linked_libs.txt

set +e
echo "[run] Rodando CFD, o output ficará armazenado em ./sw/out.txt"
if [ "$PROFILE" -eq 1 ]; then
  echo "[run] Profiling habilitado - rodando sob nsys (saida em ./sw/kernel_trace.nsys-rep)"
  nsys profile --trace=cuda --force-overwrite=true -o ./sw/kernel_trace --stats=true ./cfd >> ./sw/out.txt
else
  ./cfd >> ./sw/out.txt
fi
CFD_EXIT=$?
set -e

echo "[info] Coletando informações após o uso do código"
nvidia-smi -q > ./hw/gpu_full_query_after.txt

echo "$CFD_EXIT" > ./sw/cfd_exit_code.txt

if [ -f cuda_times.csv ]; then
  mv cuda_times.csv ./sw/
fi

mkdir -p resultados

if [ "$CFD_EXIT" -ne 0 ]; then
  echo "[warn] Erro de cfd, $CFD_EXIT - ./sw/ e ./hw/ tem informações de diagnóstico."
else
  echo "[end] CFD terminou de rodar com êxito, verificar pastas sw e hw para mais informações da run."
  echo "[end] Últimas linhas do output."
  tail -n 20 ./sw/out.txt

  if ls *.vts >/dev/null 2>&1; then
    FILE_COUNT=$(ls *.vts | wc -l)

    echo "[post] Criando arquivos vtk."
    ./build_vtks.sh tcc_lt 1 "$FILE_COUNT"

    mv *.vtk resultados/
  else
    echo "[warn] Nenhum arquivo .vts encontrado, pulando geração de vtk."
  fi
fi

if [ "$ON_SSH" -eq 1 ]; then
  REMOTE_DEST="tilt@nixos:~/p/inf/analise_desempenho/expr/$(basename $(pwd))"
  echo "[info] Tentando enviar resultados via rsync (melhor esforço)"
  rsync --verbose --progress --recursive --mkpath resultados sw hw "$REMOTE_DEST/" || \
    echo "[warn] rsync para $REMOTE_DEST falhou - rode 'make fetch' na maquina local para buscar os resultados."
fi

exit "$CFD_EXIT"


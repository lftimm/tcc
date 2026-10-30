#!/usr/bin/env bash

set -e

CFD=0
PROFILE=0
LABEL=""

while [[ $1 == -* ]]; do
  case $1 in
    -cfd)		CFD=1;;
    -profile)	PROFILE=1;;
    -label)		LABEL=$2; shift;;
    *)			echo "flag desconhecida: $1" >&2; exit 1;;
  esac
  shift
done

TIMESTAMP=$(date +%d_%m-%H:%M)
if [[ -n "$LABEL" ]]; then
  EXPR_FOLDER="expr/${LABEL}_${TIMESTAMP}"
else
  EXPR_FOLDER="expr/${TIMESTAMP}"
fi
EXPR_SRC="$EXPR_FOLDER/src"

mkdir -p $EXPR_SRC

if [[ $CFD -eq 1 ]]; then
  set -- "$@" src/build_vtks.sh src/interface.f90 src/make_vtk.f90
fi

process_file() {
  local file=$1
  local filename=$(basename "$file")
  local extension=${file##*.}

  if [[ $filename == "build_vtks.sh" ]]; then
    cp $file $EXPR_FOLDER
    return
  fi

  case $extension in
    geo)	cp $file $EXPR_FOLDER/tcc_lt.geo;;
    phys) 	cp $file $EXPR_FOLDER/tcc_lt.phys;;
    dat)	cp $file $EXPR_FOLDER/INTERFACE.dat;;
    cuf)	cp $file $EXPR_SRC/cuda_cfd.cuf;;
    *)		cp $file $EXPR_SRC;;
  esac
}

for file in "$@"; do
  process_file "$file"
done

cp Makefile run.sh $EXPR_FOLDER

if [[ $PROFILE -eq 1 ]]; then
  cp timed.slurm $EXPR_FOLDER/tupi.slurm
else
  cp tupi.slurm $EXPR_FOLDER/tupi.slurm
fi

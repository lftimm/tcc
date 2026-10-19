#!/usr/bin/bash

set -e

if [ "$#" -ne 2 ]; then
	echo "uso: $0 NOME_PROJETO PRIMEIRO_ARQUIVO NUMERO_ARQUIVOS" >&2;
	exit 1;
fi

NOME_PROJETO="$1"
PRIMEIRO_ARQUIVO="$2"
NUMERO_ARQUIVOS="$3"

for i in $(seq -f "%03g" $PRIMEIRO_ARQUIVO $NUMERO_ARQUIVOS)
do
	echo $i 
	./vtk $NOME_PROJETO $i 
done


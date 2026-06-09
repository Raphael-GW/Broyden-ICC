#!/bin/bash
GRUPOS="MEM L2CACHE FLOPS_DP FLOPS_AVX"

for grupo in $GRUPOS; do
  echo "=== Grupo: $grupo ===" >> resultados_$grupo.txt
  for N in 32 64 128 256 512 1000 2000 4000 8000 9000 10000 20000; do
    echo "--- N=$N ---" >> resultados_$grupo.txt
    likwid-perfmon -C 0 -g $grupo -m ./broyden < entradas/entrada_N${N}.txt \
      >> resultados_$grupo.txt 2>&1
  done
done
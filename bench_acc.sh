#!/bin/bash

# 1. Compilación (Necesitas el compilador nvc de NVIDIA HPC SDK)
# -acc: activa OpenACC
# -Minfo=accel: muestra información de la paralelización
# nvc -acc  -o sheet_acc sheet_acc.c -lm

# 2. Configuración de parámetros
# Formato "N|B"
PARAM_COMBIS="400|0.002 800|0.002 800|0.005 1000|0.001 1200|0.002"

echo "Iniciando Benchmark OpenACC..."
echo "N,B,Tiempo(s)" > resultados_acc.csv

for COMBI in $PARAM_COMBIS
do
    # Separar N y B
    IFS='|' read -r N_VAL B_VAL <<< "$COMBI"
    
    echo "Probando N=$N_VAL, B=$B_VAL..."

    # Medir tiempo
    T_INI=$(date +%s.%N)
    ./sheet_acc $N_VAL $B_VAL
    T_FIN=$(date +%s.%N)

    DURACION=$(echo "$T_FIN - $T_INI" | bc)
    
    # Guardar en CSV
    echo "$N_VAL,$B_VAL,$DURACION" >> resultados_acc.csv
done

echo "Hecho. Revisa resultados_acc.csv"

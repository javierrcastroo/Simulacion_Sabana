#!/bin/bash

# Combinaciones de parámetros (N y B)
# Formato: "N B"
PARAM_COMBIS="400|0.002 800|0.002 800|0.005 1000|0.001 1200|0.002"

echo "Iniciando tests..." > tiempos_secuencial.txt

for COMBI in $PARAM_COMBIS
do
    # Separar N y B usando IFS
    IFS='|' read -r N_VAL B_VAL <<< "$COMBI"
    
    echo "------------------------------------------------" >> tiempos_secuencial.txt
    echo "TEST: N=$N_VAL, B=$B_VAL" >> tiempos_secuencial.txt
    echo "------------------------------------------------" >> tiempos_secuencial.txt
    echo "Ejecutando N=$N_VAL, B=$B_VAL"
    
    # Medimos el tiempo de ejecución real
    START_TIME=$(date +%s.%N)
    ./sheet $N_VAL $B_VAL
    END_TIME=$(date +%s.%N)
    
    ELAPSED=$(echo "$END_TIME - $START_TIME" | bc)
    
    echo "Hilos: $PROC | Tiempo: $ELAPSED s" >> tiempos_secuencial.txt
done

echo "Hecho. Revisa tiempos_secuencial.txt"

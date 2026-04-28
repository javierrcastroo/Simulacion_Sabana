#!/bin/bash

# 1. Configuración de procesos y parámetros
PROCESOS="20 24 32"
PARAM_COMBIS="400|0.002 800|0.002 800|0.005 1000|0.001 1200|0.002"
EXTRA_COMBI="1500|0.001"
HOSTS="hostfile"

echo "Resultados MPI en Cluster" > resultados_cluster.txt

# Bucle principal de combinaciones normales
for COMBI in $PARAM_COMBIS
do
    IFS='|' read -r N_VAL B_VAL <<< "$COMBI"
    echo "" >> resultados_cluster.txt
    echo "--- TEST N=$N_VAL B=$B_VAL ---" >> resultados_cluster.txt
    echo "Ejecutando N=$N_VAL, B=$B_VAL"

    for NP in $PROCESOS
    do
        echo "   -> Ejecutando con $NP procesos..."
        
        START_TIME=$(date +%s.%N)
        mpirun --hostfile $HOSTS --oversubscribe -np $NP ./sheet_cluster $N_VAL $B_VAL
        END_TIME=$(date +%s.%N)
        
        ELAPSED=$(echo "$END_TIME - $START_TIME" | bc)
        echo "Procesos: $NP | Tiempo: $ELAPSED s" >> resultados_cluster.txt
    done
done

# --- SECCIÓN EXTRA ---
# Esto solo se ejecuta una vez para los procesos 24 y 32 con la combi de 1500
echo "" >> resultados_cluster.txt
echo "--- TEST EXTRA (COMBINACIÓN ESPECIAL) ---" >> resultados_cluster.txt
IFS='|' read -r NEXTRA BEXTRA <<< "$EXTRA_COMBI"

for NP in 20 24 32
do
    echo "   -> [EXTRA] Ejecutando N=$NEXTRA con $NP procesos..."
    
    START_TIME_EX=$(date +%s.%N)
    mpirun --hostfile $HOSTS --oversubscribe -np $NP ./sheet_cluster $NEXTRA $BEXTRA
    END_TIME_EX=$(date +%s.%N)
    
    ELAPSED_EX=$(echo "$END_TIME_EX - $START_TIME_EX" | bc)
    echo "Procesos: $NP (EXTRA N=$NEXTRA B=$BEXTRA) | Tiempo: $ELAPSED_EX s" >> resultados_cluster.txt
done

echo "---------------------------------------"
echo "Benchmark finalizado. Resultados en resultados_cluster.txt"

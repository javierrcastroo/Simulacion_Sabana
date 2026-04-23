#!/bin/bash

# 1. Lista de procesos a probar
# OJO: Si tienes 4 máquinas y cada una tiene, por ejemplo, 4 núcleos, 
# puedes probar hasta 16 procesos de forma eficiente.
PROCESOS="1 2 4 8 16 20"

# 2. Combinaciones N|B
PARAM_COMBIS="400|0.002 800|0.002 800|0.005 1000|0.001 1200|0.002"

# 3. Archivo de hosts (Asegúrate de que este archivo existe y tiene las IPs)
HOSTS="hostfile"

echo "Resultados MPI en Cluster" > resultados_cluster.txt

# Compilamos antes de empezar para asegurar que todos usen la última versión
# mpicc -O3 paramsEntrada_sabana_mpi_optimiced.c -o sheet_mpi -lm

for COMBI in $PARAM_COMBIS
do
    IFS='|' read -r N_VAL B_VAL <<< "$COMBI"
    echo "" >> resultados_cluster.txt
    echo "--- TEST N=$N_VAL B=$B_VAL ---" >> resultados_cluster.txt
    echo "N=$N_VAL, B=$B_VAL"

    for NP in $PROCESOS
    do
        echo "   -> Ejecutando con $NP procesos en el cluster..."
        
        # Medimos solo la ejecución
        START_TIME=$(date +%s.%N)
        
        # USAMOS --hostfile para que MPI sepa que debe ir a los otros ordenadores
        # --oversubscribe permite lanzar más procesos que núcleos físicos si fuera necesario
       mpirun --hostfile $HOSTS -x DISPLAY="" --oversubscribe -np $NP ./sheet_cluster $N_VAL $B_VAL
        
        END_TIME=$(date +%s.%N)
        
        ELAPSED=$(echo "$END_TIME - $START_TIME" | bc)
        echo "Procesos: $NP | Tiempo: $ELAPSED s" >> resultados_cluster.txt
    done
done

echo "---------------------------------------"
echo "Benchmark finalizado. Resultados en resultados_cluster.txt"

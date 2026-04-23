#!/bin/bash
PROCESOS="1 2 4 8 16 20"
# Combinaciones N|B
PARAM_COMBIS="400|0.002 800|0.002 800|0.005 1000|0.001 1200|0.002"

echo "Resultados MPI" > resultados_mpi.txt

for COMBI in $PARAM_COMBIS
do
    IFS='|' read -r N_VAL B_VAL <<< "$COMBI"
    echo "--- TEST N=$N_VAL B=$B_VAL ---" >> resultados_mpi.txt

    for NP in $PROCESOS
    do
        echo "Ejecutando N=$N_VAL con $NP procesos..."
        
        START_TIME=$(date +%s.%N)
        # Ejecución con mpirun
        mpirun --hostfile -np $NP ./sheet_mpi $N_VAL $B_VAL
        END_TIME=$(date +%s.%N)
        
        ELAPSED=$(echo "$END_TIME - $START_TIME" | bc)
        echo "Procesos: $NP | Tiempo: $ELAPSED s" >> resultados_mpi.txt
    done
done

echo "Hecho. Revisa resultados_mpi.txt"

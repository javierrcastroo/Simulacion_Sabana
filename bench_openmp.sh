#!/bin/bash

# Configuración
HILOS="1 2 4 8 16 20"
# Combinaciones de parámetros (N y B)
# Formato: "N B"
PARAM_COMBIS="400|0.002 800|0.002 800|0.005 1000|0.001 1200|0.002"

echo "Iniciando tests..." > resultados.txt

for COMBI in $PARAM_COMBIS
do
    # Separar N y B usando IFS
    IFS='|' read -r N_VAL B_VAL <<< "$COMBI"
    
    echo "------------------------------------------------" >> resultados.txt
    echo "TEST: N=$N_VAL, B=$B_VAL" >> resultados.txt
    echo "------------------------------------------------" >> resultados.txt

    for PROC in $HILOS
    do
        echo "Ejecutando N=$N_VAL con $PROC hilos..."
        export OMP_NUM_THREADS=$PROC
        
        # Medimos el tiempo de ejecución real
        # Redirigimos stderr de 'time' a un temporal
        START_TIME=$(date +%s.%N)
        ./sheet_omp $N_VAL $B_VAL
        END_TIME=$(date +%s.%N)
        
        ELAPSED=$(echo "$END_TIME - $START_TIME" | bc)
        
        echo "Hilos: $PROC | Tiempo: $ELAPSED s" >> resultados.txt
    done
done

echo "Hecho. Revisa resultados.txt"

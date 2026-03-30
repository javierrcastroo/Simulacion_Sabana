#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>

// mpicc -O2 -o sabana_mpi sabana_mpi.c -lm
// mpirun -np 4 ./sabana_mpi

#define L        10.0
#define RHO      0.1
#define T_TENS   10.0
#define B        0.002
#define AMP      0.01
#define N        1000
#define T_SIM    10.0
#define FPS      30

static double u_prev[N][N];
static double u_curr[N][N];
static double u_next[N][N];

int main(int argc, char *argv[]) {
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* 1. CONSTANTES PARA TODOS (Todos deben saber las matemáticas) */
    const double c    = sqrt(T_TENS / RHO);
    const double beta = B / RHO; 
    const double h    = L / (N - 1);
    const double dt   = 0.5 * h / c;
    const double r2   = (c * dt / h) * (c * dt / h);
    
    if (beta * dt >= 1.0) {
        if(rank == 0) printf("ERROR: damping too large\n");
        MPI_Finalize();
        return 1;
    }

    const int nstep           = (int)ceil(T_SIM / dt);
    const int steps_per_frame = (int)ceil(1.0 / (FPS * dt));
    const int nframes         = nstep / steps_per_frame + 1;
    const double omega        = (c * M_PI / L) * 2.3;
    const double c1           = 2.0 - beta * dt; 
    const double c2           = beta * dt - 1.0; 

    /* 2. INICIALIZACIÓN EXCLUSIVA DEL JEFE (rank 0) */
    FILE *fout = NULL;
    float *frame_buf = NULL;
    double *global_buf = NULL; // Buffer gigante para el MPI_Gatherv del Jefe

    if(rank == 0) {
        printf("Simulando en modo Maestro-Trabajador con %d procesos...\n", size);
        fout = fopen("frames.bin", "wb");
        int32_t hdr[2] = { (int32_t)nframes, (int32_t)N };
        fwrite(hdr, sizeof(int32_t), 2, fout);

        memset(u_prev, 0, sizeof(u_prev));
        memset(u_curr, 0, sizeof(u_curr));
        memset(u_next, 0, sizeof(u_next));

        frame_buf = malloc(N * N * sizeof(float));
        global_buf = malloc(N * N * sizeof(double));
    }

    /* 3. REPARTO DE TAREAS DINÁMICO (Soporta cualquier número de procesos) */
    int resto = N % size;
    // Si sobran filas (resto), se las repartimos a los primeros procesos
    int filas_por_proceso = (N / size) + (rank < resto ? 1 : 0);
    int j_start = rank * (N / size) + (rank < resto ? rank : resto);
    int j_end   = j_start + filas_por_proceso - 1;
    
    /* 4. PREPARAR ARRAYS PARA EL MPI_GATHERV */
    int *recvcounts = malloc(size * sizeof(int));
    int *displs = malloc(size * sizeof(int));
    int offset = 0;
    
    for (int p = 0; p < size; p++) {
        int filas_p = (N / size) + (p < resto ? 1 : 0);
        recvcounts[p] = filas_p * N; // Cuántos doubles enviará el proceso p
        displs[p] = offset;          // Dónde debe colocarlos el Jefe en el buffer
        offset += recvcounts[p];
    }

    // Buffer local para que cada proceso empaquete su trozo calculado
    double *local_buf = malloc(filas_por_proceso * N * sizeof(double));
    int frame_index = 0;

    /* BUCLE TEMPORAL */
    for (int k = 1; k <= nstep; k++) {
        double t = k * dt;

        /* A) El Jefe comparte el estado actual con todos para que puedan leer vecinos */
        MPI_Bcast(u_curr, N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(u_prev, N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        /* B) Cada trabajador procesa SOLO sus filas asignadas */
        for (int j = j_start; j <= j_end; j++) {
            if (j == 0 || j == N - 1) continue; // Saltamos los bordes superior e inferior

            for (int i = 0; i < N; i++) {
                int il = (i == 0)     ? 1     : i - 1;
                int ir = (i == N - 1) ? N - 2 : i + 1;
                
                double lap = u_curr[ir][j] + u_curr[il][j]
                           + u_curr[i][j+1] + u_curr[i][j-1]
                           - 4.0 * u_curr[i][j];
                
                u_next[i][j] = c1 * u_curr[i][j]
                             + c2 * u_prev[i][j]
                             + r2 * lap;
            }
        }

        /* Empaquetamos los datos en local_buf porque la memoria de [i][j] (columnas) no es contigua */
        int idx = 0;
        for (int j = j_start; j <= j_end; j++) {
            for (int i = 0; i < N; i++) {
                local_buf[idx++] = u_next[i][j];
            }
        }

        /* C) Recolectar las piezas con Gatherv (Soporta trozos de distinto tamaño) */
        MPI_Gatherv(local_buf, filas_por_proceso * N, MPI_DOUBLE,
                    global_buf, recvcounts, displs, MPI_DOUBLE,
                    0, MPI_COMM_WORLD);

        /* D) Tareas secuenciales del Jefe (bordes, rotar y guardar) */
        if (rank == 0) {
            
            // Desempaquetamos el puzzle gigante en u_next
            idx = 0;
            for (int j = 0; j < N; j++) {
                for (int i = 0; i < N; i++) {
                    u_next[i][j] = global_buf[idx++];
                }
            }

            // STEP 2: Borde fijo y = 0
            for (int i = 0; i < N; i++)
                u_next[i][0] = 0.0;

            // STEP 3: Borde forzado y = N-1
            for (int i = 0; i < N; i++)
                u_next[i][N-1] = AMP * sin(omega * t);

            // STEP 4: Rotar matrices
            memcpy(u_prev, u_curr, sizeof(u_curr));
            memcpy(u_curr, u_next, sizeof(u_next));

            // Guardar frame
            if (k % steps_per_frame == 0) {
                for (int j = 0; j < N; j++)
                    for (int i = 0; i < N; i++)
                        frame_buf[j * N + i] = (float)u_curr[i][j];
                
                fwrite(frame_buf, sizeof(float), N * N, fout);
                frame_index++;

                if (frame_index % 30 == 0) {
                    printf("  saved frame %d / %d\r", frame_index, nframes);
                    fflush(stdout);
                }
            }
        }
    }

    /* LIMPIEZA */
    free(recvcounts);
    free(displs);
    free(local_buf);

    if (rank == 0) {
        printf("\nDone. %d frames written to frames.bin\n", frame_index);
        free(frame_buf);
        free(global_buf);
        fclose(fout);
    }

    MPI_Finalize();
    return 0;
}
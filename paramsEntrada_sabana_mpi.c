#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>

#define L        10.0
#define RHO      0.1
#define T_TENS   10.0
#define AMP      0.01
#define T_SIM    5.0
#define FPS      30

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 3) {
        if (rank == 0) fprintf(stderr, "Uso: mpirun -np <p> %s <N> <B>\n", argv[0]);
        MPI_Finalize();
        return 1;
    }

    const int N = atoi(argv[1]);
    const double B_val = atof(argv[2]);

    const double c    = sqrt(T_TENS / RHO);
    const double beta = B_val / RHO;
    const double h    = L / (N - 1);
    const double dt   = 0.5 * h / c;
    const double r2   = (c * dt / h) * (c * dt / h);
    const double omega = (c * M_PI / L) * 2.3;

    if (beta * dt >= 1.0) {
        MPI_Finalize();
        return 1;
    }

    const int nstep = (int)ceil(T_SIM / dt);
    const double c1 = 2.0 - beta * dt;
    const double c2 = beta * dt - 1.0;

    // Reserva dinámica usando punteros a arrays para mantener sintaxis u[i][j]
    typedef double (*matrix_t)[N];
    matrix_t u_prev = calloc(1, sizeof(double[N][N]));
    matrix_t u_curr = calloc(1, sizeof(double[N][N]));
    matrix_t u_next = calloc(1, sizeof(double[N][N]));
    
    // Buffer global para recolectar datos (solo el rank 0 lo necesita realmente grande)
    double *global_buf = (rank == 0) ? malloc(N * N * sizeof(double)) : NULL;

    // Reparto de filas
    int resto = N % size;
    int filas_por_proceso = (N / size) + (rank < resto ? 1 : 0);
    int j_start = rank * (N / size) + (rank < resto ? rank : resto);
    int j_end   = j_start + filas_por_proceso - 1;

    int *recvcounts = malloc(size * sizeof(int));
    int *displs = malloc(size * sizeof(int));
    int offset = 0;
    for (int p = 0; p < size; p++) {
        int filas_p = (N / size) + (p < resto ? 1 : 0);
        recvcounts[p] = filas_p * N;
        displs[p] = offset;
        offset += recvcounts[p];
    }

    double *local_buf = malloc(filas_por_proceso * N * sizeof(double));

    for (int k = 1; k <= nstep; k++) {
        double t = k * dt;

        // A) Sincronizar estado actual (Bcast es costoso, pero mantiene tu lógica)
        MPI_Bcast(&(u_curr[0][0]), N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&(u_prev[0][0]), N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // B) Cálculo local
        for (int j = j_start; j <= j_end; j++) {
            if (j == 0 || j == N - 1) continue;
            for (int i = 0; i < N; i++) {
                int il = (i == 0) ? 1 : i - 1;
                int ir = (i == N - 1) ? N - 2 : i + 1;
                double lap = u_curr[ir][j] + u_curr[il][j] + u_curr[i][j+1] + u_curr[i][j-1] - 4.0 * u_curr[i][j];
                u_next[i][j] = c1 * u_curr[i][j] + c2 * u_prev[i][j] + r2 * lap;
            }
        }

        // Empaquetado para envío
        int idx = 0;
        for (int j = j_start; j <= j_end; j++)
            for (int i = 0; i < N; i++)
                local_buf[idx++] = u_next[i][j];

        // C) Recolectar resultados
        MPI_Gatherv(local_buf, filas_por_proceso * N, MPI_DOUBLE, global_buf, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // D) Solo el Maestro actualiza bordes y rota
        if (rank == 0) {
            idx = 0;
            for (int j = 0; j < N; j++)
                for (int i = 0; i < N; i++)
                    u_next[i][j] = global_buf[idx++];

            for (int i = 0; i < N; i++) {
                u_next[i][0] = 0.0;
                u_next[i][N-1] = AMP * sin(omega * t);
            }

            // Rotación por punteros (Swap)
            matrix_t tmp = u_prev;
            u_prev = u_curr;
            u_curr = u_next;
            u_next = tmp;
        }
    }

    free(recvcounts); free(displs); free(local_buf);
    free(u_prev); free(u_curr); free(u_next);
    if (rank == 0) free(global_buf);

    MPI_Finalize();
    return 0;
}

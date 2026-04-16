#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define L        10.0
#define RHO      0.1
#define T_TENS   10.0
#define AMP      0.01
#define T_SIM    5.0

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 3) {
        if (rank == 0) fprintf(stderr, "Uso: mpirun -np <p> %s <N> <B>\n", argv[0]);
        MPI_Finalize(); return 1;
    }

    const int N = atoi(argv[1]);
    const double B_val = atof(argv[2]);
    const double c = sqrt(T_TENS / RHO), beta = B_val / RHO, h = L / (N - 1);
    const double dt = 0.5 * h / c, r2 = (c * dt / h) * (c * dt / h);
    const double omega = (c * M_PI / L) * 2.3;
    const double c1 = 2.0 - beta * dt, c2 = beta * dt - 1.0;

    if (beta * dt >= 1.0) { MPI_Finalize(); return 1; }
    const int nstep = (int)ceil(T_SIM / dt);

    // Reparto de filas
    int resto = N % size;
    int local_rows = (N / size) + (rank < resto ? 1 : 0);
    
    // Usamos un bloque contiguo de memoria para evitar problemas de alineación
    double *data_prev = calloc((local_rows + 2) * N, sizeof(double));
    double *data_curr = calloc((local_rows + 2) * N, sizeof(double));
    double *data_next = calloc((local_rows + 2) * N, sizeof(double));

    // Punteros para usar sintaxis u[fila][columna]
    double **u_prev = malloc((local_rows + 2) * sizeof(double *));
    double **u_curr = malloc((local_rows + 2) * sizeof(double *));
    double **u_next = malloc((local_rows + 2) * sizeof(double *));

    for (int i = 0; i < local_rows + 2; i++) {
        u_prev[i] = &data_prev[i * N];
        u_curr[i] = &data_curr[i * N];
        u_next[i] = &data_next[i * N];
    }

    int up = (rank == 0) ? MPI_PROC_NULL : rank - 1;
    int down = (rank == size - 1) ? MPI_PROC_NULL : rank + 1;

    for (int k = 1; k <= nstep; k++) {
        double t = k * dt;

        // INTERCAMBIO DE HALOS (Ordenado para evitar Deadlock)
        // Los pares envían primero, los impares reciben primero
        if (rank % 2 == 0) {
            MPI_Send(u_curr[1],          N, MPI_DOUBLE, up,   0, MPI_COMM_WORLD);
            MPI_Recv(u_curr[0],          N, MPI_DOUBLE, up,   1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(u_curr[local_rows], N, MPI_DOUBLE, down, 1, MPI_COMM_WORLD);
            MPI_Recv(u_curr[local_rows+1], N, MPI_DOUBLE, down, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            MPI_Recv(u_curr[0],          N, MPI_DOUBLE, up,   1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(u_curr[1],          N, MPI_DOUBLE, up,   0, MPI_COMM_WORLD);
            MPI_Recv(u_curr[local_rows+1], N, MPI_DOUBLE, down, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(u_curr[local_rows], N, MPI_DOUBLE, down, 1, MPI_COMM_WORLD);
        }

        // CÁLCULO
        for (int j = 1; j <= local_rows; j++) {
            int global_j = rank * (N / size) + (rank < resto ? rank : resto) + (j - 1);
            if (global_j == 0 || global_j == N - 1) continue;

            for (int i = 0; i < N; i++) {
                int il = (i == 0) ? 1 : i - 1;
                int ir = (i == N - 1) ? N - 2 : i + 1;

                double lap = u_curr[j][ir] + u_curr[j][il] 
                           + u_curr[j+1][i] + u_curr[j-1][i] 
                           - 4.0 * u_curr[j][i];

                u_next[j][i] = c1 * u_curr[j][i] + c2 * u_prev[j][i] + r2 * lap;
            }
        }

        // BORDES
        int first_global_j = rank * (N / size) + (rank < resto ? rank : resto);
        int last_global_j = first_global_j + local_rows - 1;

        if (first_global_j == 0) 
            for (int i = 0; i < N; i++) u_next[1][i] = 0.0;
        if (last_global_j == N - 1) 
            for (int i = 0; i < N; i++) u_next[local_rows][i] = AMP * sin(omega * t);

        // ROTACIÓN (Swap de punteros de fila)
        double **tmp = u_prev; u_prev = u_curr; u_curr = u_next; u_next = tmp;
    }

    free(u_prev); free(u_curr); free(u_next);
    free(data_prev); free(data_curr); free(data_next);
    MPI_Finalize();
    return 0;
}

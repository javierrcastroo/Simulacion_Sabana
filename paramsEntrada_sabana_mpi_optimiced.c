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

    // Reparto de filas (Descomposición de dominio 1D)
    int resto = N % size;
    int local_rows = (N / size) + (rank < resto ? 1 : 0);
    
    // Reservamos espacio para mis filas + 2 filas extra (Halos/Fantasmas)
    // u[local_rows + 2][N]
    typedef double (*submatrix_t)[N];
    submatrix_t u_prev = calloc(local_rows + 2, sizeof(double[N]));
    submatrix_t u_curr = calloc(local_rows + 2, sizeof(double[N]));
    submatrix_t u_next = calloc(local_rows + 2, sizeof(double[N]));

    // Identificar vecinos
    int up = (rank == 0) ? MPI_PROC_NULL : rank - 1;
    int down = (rank == size - 1) ? MPI_PROC_NULL : rank + 1;

    for (int k = 1; k <= nstep; k++) {
        double t = k * dt;

        // --- INTERCAMBIO DE HALOS ---
        // Enviamos nuestra primera fila real al de arriba, recibimos su última en nuestro halo [0]
        MPI_Sendrecv(u_curr[1], N, MPI_DOUBLE, up, 0, 
                     u_curr[0], N, MPI_DOUBLE, up, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // Enviamos nuestra última fila real al de abajo, recibimos su primera en nuestro halo [last]
        MPI_Sendrecv(u_curr[local_rows], N, MPI_DOUBLE, down, 1, 
                     u_curr[local_rows + 1], N, MPI_DOUBLE, down, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // --- CÁLCULO INTERIOR ---
        for (int j = 1; j <= local_rows; j++) {
            // Determinar si esta fila 'j' es globalmente la 0 o la N-1
            int global_j = rank * (N / size) + (rank < resto ? rank : resto) + (j - 1);
            
            if (global_j == 0 || global_j == N - 1) continue;

            for (int i = 0; i < N; i++) {
                int il = (i == 0) ? 1 : i - 1;
                int ir = (i == N - 1) ? N - 2 : i + 1;

                // Accedemos a j+1 y j-1 sin miedo porque los halos tienen los datos del vecino
                double lap = u_curr[j][ir] + u_curr[j][il] 
                           + u_curr[j+1][i] + u_curr[j-1][i] 
                           - 4.0 * u_curr[j][i];

                u_next[j][i] = c1 * u_curr[j][i] + c2 * u_prev[j][i] + r2 * lap;
            }
        }

        // --- CONDICIONES DE CONTORNO (Solo los procesos de los extremos) ---
        int first_global_j = rank * (N / size) + (rank < resto ? rank : resto);
        int last_global_j = first_global_j + local_rows - 1;

        if (first_global_j == 0) {
            for (int i = 0; i < N; i++) u_next[1][i] = 0.0;
        }
        if (last_global_j == N - 1) {
            for (int i = 0; i < N; i++) u_next[local_rows][i] = AMP * sin(omega * t);
        }

        // --- ROTACIÓN DE PUNTEROS ---
        submatrix_t tmp = u_prev; u_prev = u_curr; u_curr = u_next; u_next = tmp;
    }

    free(u_prev); free(u_curr); free(u_next);
    MPI_Finalize();
    return 0;
}

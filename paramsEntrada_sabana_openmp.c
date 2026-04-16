#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

#define L       10.0
#define RHO     0.1
#define T_TENS  10.0
#define AMP     0.01
#define T_SIM   5.0   // Reducido un poco para que los tests no tarden horas
#define FPS     30

int main(int argc, char *argv[]) {
    // Valores por defecto
    int N = 500;
    double B_val = 0.002;

    if (argc > 1) N = atoi(argv[1]);
    if (argc > 2) B_val = atof(argv[2]);

    const double c = sqrt(T_TENS / RHO);
    const double beta = B_val / RHO;
    const double h = L / (N - 1);
    const double dt = 0.5 * h / c;
    const double r2 = (c * dt / h) * (c * dt / h);
    const double omega = (c * M_PI / L) * 2.3;

    if (beta * dt >= 1.0) {
        fprintf(stderr, "ERROR: damping too large\n");
        return 1;
    }

    int nstep = (int)ceil(T_SIM / dt);
    int steps_per_frame = (int)ceil(1.0 / (FPS * dt));
    
    // Reserva dinámica de matrices N x N
    double *u_prev = calloc(N * N, sizeof(double));
    double *u_curr = calloc(N * N, sizeof(double));
    double *u_next = calloc(N * N, sizeof(double));
    float *frame_buf = malloc(N * N * sizeof(float));

    const double c1 = 2.0 - beta * dt;
    const double c2 = beta * dt - 1.0;

    for (int k = 1; k <= nstep; k++) {
        double t = k * dt;

        #pragma omp parallel for schedule(static)
        for (int j = 1; j <= N - 2; j++) {
            for (int i = 0; i < N; i++) {
                int il = (i == 0) ? 1 : i - 1;
                int ir = (i == N - 1) ? N - 2 : i + 1;
                
                // Indexación 2D plana: [j*N + i]
                double lap = u_curr[ir*N + j] + u_curr[il*N + j]
                           + u_curr[i*N + (j+1)] + u_curr[i*N + (j-1)]
                           - 4.0 * u_curr[i*N + j];
                
                u_next[i*N + j] = c1 * u_curr[i*N + j]
                                + c2 * u_prev[i*N + j]
                                + r2 * lap;
            }
        }

        // Bordes (simplificado para el ejemplo)
        for (int i = 0; i < N; i++) {
            u_next[i*N + 0] = 0.0;
            u_next[i*N + (N-1)] = AMP * sin(omega * t);
        }

        memcpy(u_prev, u_curr, N * N * sizeof(double));
        memcpy(u_curr, u_next, N * N * sizeof(double));
    }

    free(u_prev); free(u_curr); free(u_next); free(frame_buf);
    return 0;
}

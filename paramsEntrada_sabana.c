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

int main(int argc, char *argv[])
{
    if (argc < 3) {
        fprintf(stderr, "Uso: %s <N> <B>\n", argv[0]);
        return 1;
    }

    int N = atoi(argv[1]);
    double B = atof(argv[2]);
    
    const double c    = sqrt(T_TENS / RHO);
    const double beta = B / RHO;
    const double h    = L / (N - 1);
    const double dt   = 0.5 * h / c;
    const double r2   = (c * dt / h) * (c * dt / h);
    const double omega = (c * M_PI / L) * 2.3;

    if (beta * dt >= 1.0) return 1;

    const int nstep = (int)ceil(T_SIM / dt);

    /* Reserva en el HEAP pero manteniendo sintaxis u[i][j].
       Esto no revienta el stack.
    */
    typedef double (*matrix_t)[N];
    matrix_t u_prev = calloc(1, sizeof(double[N][N]));
    matrix_t u_curr = calloc(1, sizeof(double[N][N]));
    matrix_t u_next = calloc(1, sizeof(double[N][N]));

    const double c1 = 2.0 - beta * dt;
    const double c2 = beta * dt - 1.0;

    for (int k = 1; k <= nstep; k++) {
        double t = k * dt;
      
        for (int j = 1; j <= N - 2; j++) {
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

        for (int i = 0; i < N; i++) {
            u_next[i][0] = 0.0;
            u_next[i][N-1] = AMP * sin(omega * t);
        }

        // En lugar de memcpy (lento), intercambiamos punteros (rápido)
        matrix_t tmp = u_prev;
        u_prev = u_curr;
        u_curr = u_next;
        u_next = tmp;
    }

    free(u_prev); 
    free(u_curr); 
    free(u_next);
    return 0;
}

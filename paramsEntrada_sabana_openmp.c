#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

#define L        10.0
#define RHO      0.1
#define T_TENS   10.0
#define AMP      0.01
#define T_SIM    5.0

int main(int argc, char *argv[])
{
    if (argc < 3) {
        fprintf(stderr, "Uso: %s <N> <B>\n", argv[0]);
        return 1;
    }

    const int N = atoi(argv[1]);
    const double B = atof(argv[2]);
    
    const double c    = sqrt(T_TENS / RHO);
    const double beta = B / RHO;
    const double h    = L / (N - 1);
    const double dt   = 0.5 * h / c;
    const double r2   = (c * dt / h) * (c * dt / h);
    const double omega = (c * M_PI / L) * 2.3;

    if (beta * dt >= 1.0) return 1;

    const int nstep = (int)ceil(T_SIM / dt);

    /* Reserva en el HEAP con punteros a arrays de tamaño N */
    double (*u_prev)[N] = calloc(1, sizeof(double[N][N]));
    double (*u_curr)[N] = calloc(1, sizeof(double[N][N]));
    double (*u_next)[N] = calloc(1, sizeof(double[N][N]));

    const double c1 = 2.0 - beta * dt;
    const double c2 = beta * dt - 1.0;

    for (int k = 1; k <= nstep; k++) {
        double t = k * dt;

        /* El paralelismo se aplica sobre las filas (j) */
        #pragma omp parallel for schedule(static)
        for (int j = 1; j <= N - 2; j++) {
            for (int i = 0; i < N; i++) {
                // Vecinos horizontales (misma fila j, columnas i-1 e i+1)
                int il = (i == 0)     ? 1      : i - 1;
                int ir = (i == N - 1) ? N - 2  : i + 1;

                // Cálculo del Laplaciano optimizado para caché (u[fila][columna])
                double lap = u_curr[j][ir] + u_curr[j][il]
                           + u_curr[j+1][i] + u_curr[j-1][i]
                           - 4.0 * u_curr[j][i];

                u_next[j][i] = c1 * u_curr[j][i]
                             + c2 * u_prev[j][i]
                             + r2 * lap;
            }
        }

        /* Bordes: Fila 0 y Fila N-1 */
        for (int i = 0; i < N; i++) {
            u_next[0][i] = 0.0;
            u_next[N-1][i] = AMP * sin(omega * t);
        }

        // Intercambio de punteros
        double (*tmp)[N] = u_prev;
        u_prev = u_curr;
        u_curr = u_next;
        u_next = tmp;
    }

    free(u_prev); 
    free(u_curr); 
    free(u_next);
    return 0;
}

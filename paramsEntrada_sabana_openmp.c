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
#define T_SIM    5.0   /* Tiempo reducido para agilizar los benchmarks */

int main(int argc, char *argv[])
{
    if (argc < 3) {
        fprintf(stderr, "Uso: %s <N> <B>\n", argv[0]);
        return 1;
    }

    // Leemos los parámetros de entrada
    const int N = atoi(argv[1]);
    const double B = atof(argv[2]);
    
    // Variables físicas
    const double c    = sqrt(T_TENS / RHO);
    const double beta = B / RHO;
    const double h    = L / (N - 1);
    const double dt   = 0.5 * h / c;
    const double r2   = (c * dt / h) * (c * dt / h);
    const double omega = (c * M_PI / L) * 2.3;

    if (beta * dt >= 1.0) {
        return 1; 
    }

    const int nstep = (int)ceil(T_SIM / dt);

    // DECLARACIÓN DE MATRICES COMO VLA (Variable Length Arrays)
    // Esto permite usar u[i][j] con el N que viene de fuera
    double u_prev[N][N];
    double u_curr[N][N];
    double u_next[N][N];

    // Inicialización
    memset(u_prev, 0, sizeof(u_prev));
    memset(u_curr, 0, sizeof(u_curr));
    memset(u_next, 0, sizeof(u_next));

    const double c1 = 2.0 - beta * dt;
    const double c2 = beta * dt - 1.0;

    for (int k = 1; k <= nstep; k++) {
        double t = k * dt;

        // STEP 1: Cálculo interior en paralelo
        #pragma omp parallel for schedule(static)
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

        // STEP 2 & 3: Bordes
        for (int i = 0; i < N; i++) {
            u_next[i][0] = 0.0;
            u_next[i][N-1] = AMP * sin(omega * t);
        }

        // STEP 4: Rotación de datos
        memcpy(u_prev, u_curr, sizeof(u_curr));
        memcpy(u_curr, u_next, sizeof(u_next));
    }

    return 0;
}

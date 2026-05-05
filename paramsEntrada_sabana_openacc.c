#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <openacc.h>

#define L        10.0
#define RHO      0.1
#define T_TENS   10.0
#define AMP      0.01
#define T_SIM    5.0

int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Uso: %s <N> <B>\n", argv[0]);
        return 1;
    }

    const int N = atoi(argv[1]);
    const double B_val = atof(argv[2]);

    const double c = sqrt(T_TENS / RHO), beta = B_val / RHO, h = L / (N - 1);
    const double dt = 0.5 * h / c, r2 = (c * dt / h) * (c * dt / h);
    const double omega = (c * M_PI / L) * 2.3;
    const double c1 = 2.0 - beta * dt, c2 = beta * dt - 1.0;

    if (beta * dt >= 1.0) return 1;
    const int nstep = (int)ceil(T_SIM / dt);

    typedef double (*matrix_t)[N];
    matrix_t u_prev = calloc(N, sizeof(double[N]));
    matrix_t u_curr = calloc(N, sizeof(double[N]));
    matrix_t u_next = calloc(N, sizeof(double[N]));

    // Usamos 'copy' para asegurar que los punteros iniciales estén vinculados
    #pragma acc data copy(u_prev[0:N][0:N], u_curr[0:N][0:N], u_next[0:N][0:N])
    {
        for (int k = 1; k <= nstep; k++) {
            double t = k * dt;

            // STEP 1: Cálculo interior paralelo en GPU
            #pragma acc parallel loop collapse(2) present(u_prev, u_curr, u_next)
            for (int j = 1; j <= N - 2; j++) {
                for (int i = 0; i < N; i++) {
                    int il = (i == 0)     ? 1     : i - 1;
                    int ir = (i == N - 1) ? N - 2 : i + 1;
                   
                    double lap = u_curr[j][ir] + u_curr[j][il]
                               + u_curr[j+1][i] + u_curr[j-1][i]
                               - 4.0 * u_curr[j][i];
                   
                    u_next[j][i] = c1 * u_curr[j][i]
                                 + c2 * u_prev[j][i]
                                 + r2 * lap;
                }
            }

            // STEP 2: Bordes en GPU
            #pragma acc parallel loop present(u_next)
            for (int i = 0; i < N; i++) {
                u_next[0][i] = 0.0;
                u_next[N-1][i] = AMP * sin(omega * t);
            }

            // STEP 3: Swap de punteros (Muy rápido, solo cambia la dirección)
            // Al estar dentro de 'acc data', OpenACC rastrea qué puntero host
            // corresponde a qué zona de la GPU.
            matrix_t tmp = u_prev;
            u_prev = u_curr;
            u_curr = u_next;
            u_next = tmp;
        }
       
        // CRUCIAL: Antes de salir, pedimos a la GPU que actualice el 'u_curr' actual en el Host
        #pragma acc update host(u_curr[0:N][0:N])
    }

    // Validación con formato científico para ver precisión
    int centro = N / 2;
    printf("Check: u[N-1][%d] = %e | u[%d][%d] = %e\n",
            centro, u_curr[N-1][centro], centro, centro, u_curr[centro][centro]);
   
    free(u_prev); free(u_curr); free(u_next);
    return 0;
}

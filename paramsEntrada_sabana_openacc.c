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

    // Reserva de memoria (igual que antes)
    typedef double (*matrix_t)[N];
    matrix_t u_prev = calloc(N, sizeof(double[N]));
    matrix_t u_curr = calloc(N, sizeof(double[N]));
    matrix_t u_next = calloc(N, sizeof(double[N]));

    // --- REGIÓN DE DATOS OPENACC ---
    // 'copyin': sube u_prev y u_curr a la GPU al empezar.
    // 'create': reserva espacio para u_next en la GPU sin subir basura.
    // 'copyout': baja el resultado final u_curr a la RAM al terminar.
    #pragma acc data copyin(u_prev[0:N][0:N], u_curr[0:N][0:N]) create(u_next[0:N][0:N]) copyout(u_curr[0:N][0:N])
    {
        for (int k = 1; k <= nstep; k++) {
            double t = k * dt;

            // STEP 1: Cálculo interior paralelo en GPU
            // 'collapse(2)' hace que la GPU trate los dos bucles como uno solo de N*N hilos
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

            // STEP 2 & 3: Bordes (también en la GPU para no mover datos)
            #pragma acc parallel loop present(u_next)
            for (int i = 0; i < N; i++) {
                u_next[0][i] = 0.0;
                u_next[N-1][i] = AMP * sin(omega * t);
            }

            // STEP 4: Rotación de punteros (Swap)
            // IMPORTANTE: En OpenACC, si swapeas en el host, 
            for (int j = 0; j < N; j++) {
                for (int i = 0; i < N; i++) {
                    u_prev[j][i] = u_curr[j][i]; // El presente pasa a ser pasado
                    u_curr[j][i] = u_next[j][i]; // El futuro pasa a ser el presente
                }
            }
        }
    }

    // Imprimimos un valor periférico y uno central como checksum
    printf("Check: u[N-1][N/2] = %f | u[N/2][N/2] = %f\n", u_curr[N-1][N/2], u_curr[N/2][N/2]);
    free(u_prev); free(u_curr); free(u_next);
    return 0;
}

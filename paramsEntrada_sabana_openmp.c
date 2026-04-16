#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>

#define L        10.0      /* side length [m]                          */
#define RHO      0.1      /* surface mass density [kg/m²]             */
#define T_TENS   10.0     /* surface tension [N/m]                    */
                          /* increase to damp faster, decrease for    */
                          /* longer oscillations. typical range: 0.1-2*/
#define AMP      0.01     /* forcing amplitude [m]                    */

/* ------------------------------------------------------------------ */
/*  Numerical parameters                                                */
/* ------------------------------------------------------------------ */
#define T_SIM    10.0      /* total simulation time [s]                */

/* ------------------------------------------------------------------ */
/*  Output parameters                                                   */
/* ------------------------------------------------------------------ */
#define FPS      30       /* frames per second to save                */

/* ------------------------------------------------------------------ */

static double u_prev[N][N];
static double u_curr[N][N];
static double u_next[N][N];

int main(int argc, char *argv[])
{
    if (argc < 3) {
        fprintf(stderr, "Uso: %s <N> <B>\n", argv[0]);
        return 1;
    }

    int N = atoi(argv[1]);
    double B = atof(argv[2]);
    
    const double c    = sqrt(T_TENS / RHO);
    const double beta = B / RHO;             /* damping rate [1/s]    */
    const double h    = L / (N - 1);
    const double dt   = 0.5 * h / c;
    const double r2   = (c * dt / h) * (c * dt / h);

    /* Damping stability: beta*dt must be < 1 */
    if (beta * dt >= 1.0) {
        fprintf(stderr, "ERROR: damping too large (beta*dt=%.4f >= 1)\n",
                beta * dt);
        fprintf(stderr, "       Reduce B or increase N.\n");
        return 1;
    }

    const int nstep           = (int)ceil(T_SIM / dt);
    const int steps_per_frame = (int)ceil(1.0 / (FPS * dt));
    const int nframes         = nstep / steps_per_frame + 1;

    /* Forcing frequency: mode (m=0, n=1) — fundamental */
    const double omega = (c * M_PI / L) * 2.3;

    memset(u_prev, 0, sizeof(u_prev));
    memset(u_curr, 0, sizeof(u_curr));
    memset(u_next, 0, sizeof(u_next));

    /* Precompute damping coefficients */
    const double c1 = 2.0 - beta * dt;   /* coefficient of u_curr    */
    const double c2 = beta * dt - 1.0;   /* coefficient of u_prev    */

    int frame_index = 0;

    for (int k = 1; k <= nstep; k++) {

        double t = k * dt;

        /* ---------------------------------------------------------- */
        /* STEP 1: update interior nodes                               */
        /* Each j-row is independent: safe to parallelize with OpenMP  */
        /* ---------------------------------------------------------- */
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

        /* STEP 2: fixed boundary y = 0 */
        for (int i = 0; i < N; i++)
            u_next[i][0] = 0.0;

        /* STEP 3: forced boundary y = N-1 */
        for (int i = 0; i < N; i++)
            u_next[i][N-1] = AMP * sin(omega * t);

        /* STEP 4: rotate arrays */
        memcpy(u_prev, u_curr, sizeof(u_curr));
        memcpy(u_curr, u_next, sizeof(u_next));
    }
    return 0;
}

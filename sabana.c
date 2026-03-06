/*
 * sheet.c
 * Simulation of normal vibration modes of an elastic sheet.
 *
 * Model:  d²u/dt² = c²∇²u + α∇²(du/dt) - β(du/dt) - g
 *
 * Boundary conditions:
 *   y = 0      : fixed (Dirichlet homogeneous)
 *   y = N-1    : forced A*sin(omega*t) (Dirichlet non-homogeneous)
 *   x = 0, N-1 : free (Neumann, ghost nodes)
 *
 * Output: frames.bin
 *   header: int32 nframes, int32 N
 *   then nframes blocks of N*N float32 values, indexed [row=y, col=x]
 *
 * Compile:
 *   gcc -O2 -o sheet sheet.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>

/* ------------------------------------------------------------------ */
/*  Physical parameters                                                 */
/* ------------------------------------------------------------------ */
#define L        10.0      /* side length [m]                          */
#define RHO      0.1      /* surface mass density [kg/m²]             */
#define T_TENS   10.0     /* surface tension [N/m]                    */
#define ETA      0.0002    /* internal viscosity Kelvin-Voigt [N·s/m]  */
#define B        0.1      /* external viscous damping [kg/(m²·s)]     */
#define AMP      0.01     /* forcing amplitude [m]                    */
#define G        9.81     /* gravity [m/s²]                           */

/* ------------------------------------------------------------------ */
/*  Numerical parameters                                                */
/* ------------------------------------------------------------------ */
#define N        500     /* grid nodes per side                      */
#define T_SIM    20.0      /* total simulation time [s]                */

/* ------------------------------------------------------------------ */
/*  Output parameters                                                   */
/* ------------------------------------------------------------------ */
#define FPS      30       /* frames per second to save                */

/* ------------------------------------------------------------------ */

static double u_prev[N][N];
static double u_curr[N][N];
static double u_next[N][N];
static double L_prev[N][N];
static double L_curr[N][N];

int main(void)
{
    /* Derived physical constants */
    const double c     = sqrt(T_TENS / RHO);
    const double alpha = ETA / RHO;
    const double beta  = B   / RHO;

    /* Spatial and temporal steps */
    const double h  = L / (N - 1);
    const double dt = 0.5 * h / c;

    /* Stability checks */
    const double r2  = (c * dt / h) * (c * dt / h);
    const double gam = alpha * dt / (h * h);
    if (gam > 0.25) {
        fprintf(stderr, "ERROR: parabolic CFL violated (gam=%.4f > 0.25)\n", gam);
        return 1;
    }

    const int nstep           = (int)ceil(T_SIM / dt);
    const int steps_per_frame = (int)ceil(1.0 / (FPS * dt));
    const int nframes         = nstep / steps_per_frame + 1;

    /* Forcing frequency: mode (0,1) - fundamental */
    const double omega = (c * M_PI / L) * sqrt(0.0*0.0 + 1.0*1.0);

    printf("c      = %.4f m/s\n", c);
    printf("h      = %.6f m\n",   h);
    printf("dt     = %.6f s\n",   dt);
    printf("r      = %.4f  (CFL limit: %.4f)\n", c*dt/h, 1.0/sqrt(2.0));
    printf("gamma  = %.6f  (limit: 0.25)\n", gam);
    printf("omega  = %.4f rad/s  (f = %.4f Hz)\n", omega, omega/(2*M_PI));
    printf("nstep  = %d\n",   nstep);
    printf("nframes= %d\n",   nframes);

    /* Open output binary file */
    FILE *fout = fopen("frames.bin", "wb");
    if (!fout) { perror("frames.bin"); return 1; }

    /* Header: nframes and N as int32 */
    int32_t hdr[2] = { (int32_t)nframes, (int32_t)N };
    fwrite(hdr, sizeof(int32_t), 2, fout);

    /* Initialize arrays */
    memset(u_prev, 0, sizeof(u_prev));
    memset(u_curr, 0, sizeof(u_curr));
    memset(u_next, 0, sizeof(u_next));
    memset(L_prev, 0, sizeof(L_prev));
    memset(L_curr, 0, sizeof(L_curr));

    /* Buffer to convert one frame to float32 for writing */
    float *frame_buf = malloc(N * N * sizeof(float));
    if (!frame_buf) { fprintf(stderr, "malloc failed\n"); return 1; }

    int frame_index = 0;

    for (int k = 1; k <= nstep; k++) {

        double t = k * dt;

        /* STEP 1: laplacian numerator with Neumann BC */
        for (int j = 1; j <= N - 2; j++) {
            for (int i = 0; i < N; i++) {
                int il = (i == 0)     ? 1     : i - 1;
                int ir = (i == N - 1) ? N - 2 : i + 1;
                L_curr[i][j] = u_curr[ir][j] + u_curr[il][j]
                             + u_curr[i][j+1] + u_curr[i][j-1]
                             - 4.0 * u_curr[i][j];
            }
        }

        /* STEP 2: update interior nodes */
        for (int j = 1; j <= N - 2; j++) {
            for (int i = 0; i < N; i++) {
                u_next[i][j] = (2.0 - beta * dt) * u_curr[i][j]
                             + (beta * dt - 1.0)  * u_prev[i][j]
                             + (r2 + gam)          * L_curr[i][j]
                             -  gam                * L_prev[i][j]
                             -  G * dt * dt;
            }
        }

        /* STEP 3: fixed boundary y = 0 */
        for (int i = 0; i < N; i++)
            u_next[i][0] = 0.0;

        /* STEP 4: forced boundary y = N-1 */
        for (int i = 0; i < N; i++)
            u_next[i][N-1] = AMP * sin(omega * t);

        /* STEP 5: rotate arrays */
        memcpy(L_prev, L_curr, sizeof(L_curr));
        memcpy(u_prev, u_curr, sizeof(u_curr));
        memcpy(u_curr, u_next, sizeof(u_next));

        /* Save frame */
        if (k % steps_per_frame == 0) {
            /* store as frame_buf[j*N + i] so Python reads u[y, x] */
            for (int j = 0; j < N; j++)
                for (int i = 0; i < N; i++)
                    frame_buf[j * N + i] = (float)u_curr[i][j];
            fwrite(frame_buf, sizeof(float), N * N, fout);
            frame_index++;

            if (frame_index % 30 == 0)
                printf("  saved frame %d / %d\r", frame_index, nframes);
            fflush(stdout);
        }
    }

    printf("\nDone. %d frames written to frames.bin\n", frame_index);

    free(frame_buf);
    fclose(fout);
    return 0;
}

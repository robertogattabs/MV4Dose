#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include <complex.h>
#include <limits.h>
#define min(a, b) (((a) < (b)) ? (a) : (b));



/* Ritorna 1 se axis è strettamente crescente, altrimenti 0 */
static int is_strictly_increasing(const double *axis, int n) {
    for (int i = 1; i < n; ++i) if (!(axis[i] > axis[i-1])) return 0;
    return 1;
}

/* Ricerca binaria: trova i tale che axis[i] <= v < axis[i+1]
 Ritorna -1 se v < axis[0] o v > axis[n-1].
 Caso v == axis[n-1]: ritorna n-2 e imposta *t = 1.0 */
static int bracket_index(double v, const double *axis, int n, double *t) {
    if (v < axis[0] || v > axis[n-1]) return -1;
    if (v == axis[n-1]) {                 /* esattamente al bordo alto */
*t = 1.0;
        return n - 2;
    }
    int lo = 0, hi = n - 2;               /* ultima cella è n-2 .. n-1 */
while (lo <= hi) {
    int mid = (lo + hi) >> 1;
    double a0 = axis[mid], a1 = axis[mid + 1];
    if (v < a0) {
        hi = mid - 1;
    } else if (v > a1) {
        lo = mid + 1;
    } else {
        *t = (a1 != a0) ? (v - a0) / (a1 - a0) : 0.0;
        return mid;
    }
}
return -1; /* non dovrebbe capitare con assi crescenti */
}

/* newXps è linearizzato z-veloce: idx = i*ny*nz + j*nz + k */
void interpola_suAltraMatricePunti(
        /* griglia destinazione (CT) */
        double *xPos_CT, double *yPos_CT, double *zPos_CT,
        /* griglia sorgente (DS) */
        double *xPos_DS, double *yPos_DS, double *zPos_DS,
        /* lunghezze */
        int *len_x_CT, int *len_y_CT, int *len_z_CT,
        int *len_x_DS, int *len_y_DS, int *len_z_DS,
        /* valori sorgente linearizzati (z-veloce) e output */
        double *newXps,
        double *returnMatrix)
{
    const int nxCT = *len_x_CT, nyCT = *len_y_CT, nzCT = *len_z_CT;
    const int nxDS = *len_x_DS, nyDS = *len_y_DS, nzDS = *len_z_DS;

    /* controlli minimi */
    if (nxDS < 2 || nyDS < 2 || nzDS < 2) {
        Rprintf("[interp] Griglia DS troppo piccola (>=2 per asse richiesto).\n");
        R_FlushConsole();
        /* output già 0 se allocato con double(...); in ogni caso azzero */
        const size_t N = (size_t)nxCT * nyCT * nzCT;
        for (size_t q = 0; q < N; ++q) returnMatrix[q] = 0.0;
        return;
    }
    if (!is_strictly_increasing(xPos_DS, nxDS) ||
        !is_strictly_increasing(yPos_DS, nyDS) ||
        !is_strictly_increasing(zPos_DS, nzDS)) {
        Rprintf("[interp] Assi DS non strettamente crescenti. Restituisco 0.\n");
        R_FlushConsole();
        const size_t N = (size_t)nxCT * nyCT * nzCT;
        for (size_t q = 0; q < N; ++q) returnMatrix[q] = 0.0;
        return;
    }

    size_t out_idx = 0;
    for (int ix = 0; ix < nxCT; ++ix) {
        const double x = xPos_CT[ix];
        double tx; int i = bracket_index(x, xPos_DS, nxDS, &tx);

        for (int iy = 0; iy < nyCT; ++iy) {
            const double y = yPos_CT[iy];
            double ty; int j = bracket_index(y, yPos_DS, nyDS, &ty);

            for (int iz = 0; iz < nzCT; ++iz) {
                const double z = zPos_CT[iz];
                double tz; int k = bracket_index(z, zPos_DS, nzDS, &tz);

                if (i < 0 || j < 0 || k < 0) {
                    /* fuori dal dominio DS: set a 0 */
                    returnMatrix[out_idx++] = 0.0;
                    continue;
                }

                /* indici vertici nel volume DS (z-veloce) */
                const int strideYZ = nyDS * nzDS;
                const int strideZ  = nzDS;

                const int idx000 =  i      * strideYZ +  j      * strideZ +  k;
                const int idx001 =  i      * strideYZ +  j      * strideZ + (k+1);
                const int idx010 =  i      * strideYZ + (j+1)   * strideZ +  k;
                const int idx011 =  i      * strideYZ + (j+1)   * strideZ + (k+1);
                const int idx100 = (i + 1) * strideYZ +  j      * strideZ +  k;
                const int idx101 = (i + 1) * strideYZ +  j      * strideZ + (k+1);
                const int idx110 = (i + 1) * strideYZ + (j+1)   * strideZ +  k;
                const int idx111 = (i + 1) * strideYZ + (j+1)   * strideZ + (k+1);

                const double v000 = newXps[idx000];
                const double v001 = newXps[idx001];
                const double v010 = newXps[idx010];
                const double v011 = newXps[idx011];
                const double v100 = newXps[idx100];
                const double v101 = newXps[idx101];
                const double v110 = newXps[idx110];
                const double v111 = newXps[idx111];

                /* trilineare: interpolo lungo x, poi y, poi z */
                const double c00 = v000*(1.0 - tx) + v100*tx;
                const double c01 = v001*(1.0 - tx) + v101*tx;
                const double c10 = v010*(1.0 - tx) + v110*tx;
                const double c11 = v011*(1.0 - tx) + v111*tx;

                const double c0  = c00*(1.0 - ty) + c10*ty;
                const double c1  = c01*(1.0 - ty) + c11*ty;

                const double c   = c0*(1.0 - tz) + c1*tz;

                returnMatrix[out_idx++] = c;
            }
        }
    }
    /* opzionale: messaggio di fine */
    /* Rprintf("[interp] Completata.\n"); R_FlushConsole(); */
}

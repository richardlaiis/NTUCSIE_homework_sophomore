#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include "cblas.h"

void steepest_descent(int d, double *A, double *b, double *x, double tol) {
    double *r = (double *)malloc(d * sizeof(double));
    double *q = (double *)malloc(d * sizeof(double));
    
    // x = 0
    for (int i = 0; i < d; i++) {
        x[i] = 0.0;
    }
    
    // r = b
    cblas_dcopy(d, b, 1, r, 1);
    
    double eps = cblas_dnrm2(d, r, 1);
    
    while (eps > tol) {
        // q = A * r
        cblas_dgemv(CblasRowMajor, CblasNoTrans, d, d, 1.0, A, d, r, 1, 0.0, q, 1);
        
        // alpha = (r^T r) / (r^T q)
        double rTr = cblas_ddot(d, r, 1, r, 1);
        double rTq = cblas_ddot(d, r, 1, q, 1);
        
        // Prevent division by zero if r is already exactly 0 but loop continued
        if (rTq == 0) break;
        double alpha = rTr / rTq;
        
        // x = x + alpha * r
        cblas_daxpy(d, alpha, r, 1, x, 1);
        
        // r = b - A * x
        cblas_dcopy(d, b, 1, r, 1);
        cblas_dgemv(CblasRowMajor, CblasNoTrans, d, d, -1.0, A, d, x, 1, 1.0, r, 1);
        
        // eps = ||r||_2
        eps = cblas_dnrm2(d, r, 1);
    }
    
    free(r);
    free(q);
}

// Generate SPD matrix
// A_ij = 1.0 if i != j, A_ii = d + 1.0
// b_i = i + 1.0
void generate_spd_matrix(int d, double **A, double **b) {
    *A = (double *)malloc(d * d * sizeof(double));
    *b = (double *)malloc(d * sizeof(double));
    
    for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) {
            if (i == j) {
                (*A)[i * d + j] = d + 1.0;
            } else {
                (*A)[i * d + j] = 1.0;
            }
        }
        (*b)[i] = i + 1.0;
    }
}

int main(int argc, char **argv) {
    int d = 100; // Default size
    if (argc > 1) {
        d = atoi(argv[1]);
    }

    double tol = 1e-6;
    double *A, *b;

    generate_spd_matrix(d, &A, &b);

    double *x = (double *)malloc(d * sizeof(double));

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    steepest_descent(d, A, b, x, tol);

    clock_gettime(CLOCK_MONOTONIC, &end);
    double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("Dimension: %d x %d\n", d, d);
    printf("Elapsed time: %f seconds\n", elapsed);

    printf("Solution for first 5 elements:\n");
    for (int i = 0; i < 5 && i < d; i++) {
        printf("x[%d] = %f\n", i, x[i]);
    }

    free(A);
    free(b);
    free(x);

    return 0;
}
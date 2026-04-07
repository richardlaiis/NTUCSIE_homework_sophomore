// my_steepest_descent.c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "cblas.h"
#define tol 1e-8
#define seed 53

double random_0_to_1() {
    return (double)rand() / (double)RAND_MAX;
}

void generate_spd(double *A, int n) {
    double *M = (double *)malloc((n*n) * sizeof(double));
    for (int i = 0; i < n*n; i++) M[i] = random_0_to_1();
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, n, n, 1.0, M, n, M, n, 0.0, A, n);
    for (int i = 0; i < n; i++) A[i*n+i] += (double)n;
    free(M);
}

void generate_vector(double *v, int n) {
    for (int i = 0; i < n; i++) v[i] = random_0_to_1();
}

void steepest_descent(double *A, double *b, double *x, int n) {
    double *r = (double *)malloc(n * sizeof(double));
    double *Ar = (double *)malloc(n * sizeof(double));
    double *Ax = (double *)malloc(n * sizeof(double));

    memset(x, 0, n*sizeof(double));
    cblas_dcopy(n, b, 1, r, 1);
    double eps = cblas_dnrm2(n, r, 1);
    while (eps > tol) {
        cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0, A, n, r, 1, 0.0, Ar, 1);
        double rr = cblas_ddot(n, r, 1, r, 1);
        double rAr = cblas_ddot(n, r, 1, Ar, 1);
        double alpha = rr / rAr;
        cblas_daxpy(n, alpha, r, 1, x, 1);
        cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0, A, n, x, 1, 0.0, Ax, 1);
        cblas_dcopy(n, b, 1, r, 1);
        cblas_daxpy(n, -1.0, Ax, 1, r, 1);
        eps = cblas_dnrm2(n, r, 1);
    }

    free(r);
    free(Ar);
    free(Ax);
}

void write_vector(char *path, double *v, int n) {
    FILE *fp = fopen(path, "w");
    for (int i = 0; i < n; i++) {
        fprintf(fp, "%f\n", v[i]);
    }
    fclose(fp);
}

void write_matrix(char *path, double *A, int n) {
    FILE *fp = fopen(path, "w");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fprintf(fp, "%f%c", A[i*n+j], " \n"[j==n-1]);
        }
    }
    fclose(fp);
}

int main(int argc, char **argv) {
    int n = 100;
    if (argc > 1) {
        n = atoi(argv[1]);
    }
    srand(seed);

    double *A = (double *)malloc((n*n) * sizeof(double));
    double *b = (double *)malloc(n * sizeof(double));
    double *x = (double *)malloc(n * sizeof(double));

    generate_spd(A, n);
    generate_vector(b, n);
    steepest_descent(A, b, x, n);

    double *Ax = (double *)malloc(n * sizeof(double));
    cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0, A, n, x, 1, 0.0, Ax, 1);

    write_matrix("A", A, n);
    write_vector("x", x, n);
    write_vector("b", b, n);


}
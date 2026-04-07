#include <cblas.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct {
	int iterations;
	double residual_norm;
	int converged;
} SDResult;

static double uniform_random_0_1(void) {
	return (double)rand() / (double)RAND_MAX;
}

static void generate_spd_matrix(double *A, int n) {
	double *M = (double *)malloc((size_t)n * (size_t)n * sizeof(double));
	if (M == NULL) {
		fprintf(stderr, "Failed to allocate memory for temporary matrix M.\n");
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < n * n; ++i) {
		M[i] = uniform_random_0_1();
	}

	cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, n, n, 1.0, M, n, M,
				n, 0.0, A, n);

	for (int i = 0; i < n; ++i) {
		A[i * n + i] += (double)n;
	}

	free(M);
}

static void generate_vector(double *v, int n) {
	for (int i = 0; i < n; ++i) {
		v[i] = uniform_random_0_1();
	}
}

static SDResult steepest_descent(const double *A, const double *b, double *x, int n,
								 double tol, int max_iter) {
	double *r = (double *)malloc((size_t)n * sizeof(double));
	double *Ar = (double *)malloc((size_t)n * sizeof(double));
	double *Ax = (double *)malloc((size_t)n * sizeof(double));

	if (r == NULL || Ar == NULL || Ax == NULL) {
		fprintf(stderr, "Failed to allocate memory for steepest descent buffers.\n");
		free(r);
		free(Ar);
		free(Ax);
		exit(EXIT_FAILURE);
	}

	memset(x, 0, (size_t)n * sizeof(double));
	cblas_dcopy(n, b, 1, r, 1);

	double residual = cblas_dnrm2(n, r, 1);
	int iter = 0;

	while (residual > tol && iter < max_iter) {
		cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0, A, n, r, 1, 0.0, Ar, 1);

		const double rr = cblas_ddot(n, r, 1, r, 1);
		const double rAr = cblas_ddot(n, r, 1, Ar, 1);

		if (fabs(rAr) < 1e-20) {
			break;
		}

		const double alpha = rr / rAr;

		cblas_daxpy(n, alpha, r, 1, x, 1);

		cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0, A, n, x, 1, 0.0, Ax, 1);
		cblas_dcopy(n, b, 1, r, 1);
		cblas_daxpy(n, -1.0, Ax, 1, r, 1);

		residual = cblas_dnrm2(n, r, 1);
		++iter;
	}

	free(r);
	free(Ar);
	free(Ax);

	SDResult result;
	result.iterations = iter;
	result.residual_norm = residual;
	result.converged = (residual <= tol);
	return result;
}

static double cpu_seconds_since(clock_t start, clock_t end) {
	return (double)(end - start) / (double)CLOCKS_PER_SEC;
}

static void write_vector(const char *path, const double *v, int n) {
	FILE *fp = fopen(path, "w");
	if (fp == NULL) {
		fprintf(stderr, "Failed to open %s for writing.\n", path);
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < n; ++i) {
		fprintf(fp, "%.17g\n", v[i]);
	}

	fclose(fp);
}

static void write_matrix(const char *path, const double *A, int n) {
	FILE *fp = fopen(path, "w");
	if (fp == NULL) {
		fprintf(stderr, "Failed to open %s for writing.\n", path);
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			fprintf(fp, j + 1 == n ? "%.17g" : "%.17g ", A[i * n + j]);
		}
		fputc('\n', fp);
	}

	fclose(fp);
}

int main(int argc, char **argv) {
	int n = 100;
	double tol = 1e-8;
	int max_iter = 200000;
	unsigned int seed = 42;
	const char *dump_prefix = NULL;

	if (argc > 1) {
		n = atoi(argv[1]);
	}
	if (argc > 2) {
		tol = atof(argv[2]);
	}
	if (argc > 3) {
		max_iter = atoi(argv[3]);
	}
	if (argc > 4) {
		seed = (unsigned int)strtoul(argv[4], NULL, 10);
	}
	if (argc > 5) {
		dump_prefix = argv[5];
	}

	if (n <= 0 || tol <= 0.0 || max_iter <= 0) {
		fprintf(stderr,
				"Usage: %s [n>0] [tol>0] [max_iter>0] [seed] [dump_prefix]\n"
				"Example: %s 1000 1e-8 200000 42 verify\n",
				argv[0], argv[0]);
		return EXIT_FAILURE;
	}

	srand(seed);

	double *A = (double *)malloc((size_t)n * (size_t)n * sizeof(double));
	double *b = (double *)malloc((size_t)n * sizeof(double));
	double *x = (double *)malloc((size_t)n * sizeof(double));
	double *x_true = (double *)malloc((size_t)n * sizeof(double));

	if (A == NULL || b == NULL || x == NULL || x_true == NULL) {
		fprintf(stderr, "Failed to allocate matrix/vector memory.\n");
		free(A);
		free(b);
		free(x);
		free(x_true);
		return EXIT_FAILURE;
	}

	generate_spd_matrix(A, n);
	generate_vector(x_true, n);

	cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0, A, n, x_true, 1, 0.0, b, 1);

	clock_t t0 = clock();
	SDResult result = steepest_descent(A, b, x, n, tol, max_iter);
	clock_t t1 = clock();

	double *err = (double *)malloc((size_t)n * sizeof(double));
	if (err == NULL) {
		fprintf(stderr, "Failed to allocate error buffer.\n");
		free(A);
		free(b);
		free(x);
		free(x_true);
		return EXIT_FAILURE;
	}

	cblas_dcopy(n, x, 1, err, 1);
	cblas_daxpy(n, -1.0, x_true, 1, err, 1);
	const double rel_error = cblas_dnrm2(n, err, 1) / cblas_dnrm2(n, x_true, 1);

	printf("Steepest Descent (CBLAS)\n");
	printf("n            : %d\n", n);
	printf("tolerance    : %.3e\n", tol);
	printf("max_iter     : %d\n", max_iter);
	printf("iterations   : %d\n", result.iterations);
	printf("converged    : %s\n", result.converged ? "yes" : "no");
	printf("||r||_2      : %.6e\n", result.residual_norm);
	printf("relative err : %.6e\n", rel_error);
	printf("elapsed (s)  : %.6f\n", cpu_seconds_since(t0, t1));

	const int print_count = (n < 10) ? n : 10;
	printf("solution x (first %d entries):\n", print_count);
	for (int i = 0; i < print_count; ++i) {
		printf("x[%d] = %.10f\n", i, x[i]);
	}

	if (dump_prefix != NULL) {
		char path[512];

		snprintf(path, sizeof(path), "%s_A.txt", dump_prefix);
		write_matrix(path, A, n);

		snprintf(path, sizeof(path), "%s_b.txt", dump_prefix);
		write_vector(path, b, n);

		snprintf(path, sizeof(path), "%s_x.txt", dump_prefix);
		write_vector(path, x, n);

		snprintf(path, sizeof(path), "%s_x_true.txt", dump_prefix);
		write_vector(path, x_true, n);

		printf("dumped verification files with prefix '%s'\n", dump_prefix);
	}

	free(A);
	free(b);
	free(x);
	free(x_true);
	free(err);
	return EXIT_SUCCESS;
}

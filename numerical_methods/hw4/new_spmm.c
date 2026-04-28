#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <mkl.h>

#define ABS(x) ((x) < 0.0f ? -(x) : (x))
#define eps 1e-2

// Unified structure for both CSR and CSC formats
typedef struct { int r, c, nnz, is_csr; float *v; int *idx, *ptr; } SpMat;

// Unified SpMM function for both CSR and CSC
void spmm(SpMat *A, float *B, float *C, int n, int use_omp) {
    for (int i = 0; i < A->r * n; i++) C[i] = 0.0f;
    
    if (A->is_csr) {
        #pragma omp parallel for if(use_omp)
        for (int i = 0; i < A->r; i++) {
            for (int p = A->ptr[i]; p < A->ptr[i+1]; p++) {
                for (int j = 0; j < n; j++) 
                    C[i*n + j] += A->v[p] * B[A->idx[p]*n + j];
            }
        }
    } else {
        for (int k = 0; k < A->c; k++) {
            for (int p = A->ptr[k]; p < A->ptr[k+1]; p++) {
                for (int j = 0; j < n; j++) 
                    C[A->idx[p]*n + j] += A->v[p] * B[k*n + j];
            }
        }
    }
}

// Unified random matrix generator (NNZ = 32 * m constraint)
void gen_sp(SpMat *A, int m, int k, int is_csr) {
    int nnz = 32 * m, out_dim = is_csr ? m : k;
    A->r = m; A->c = k; A->nnz = nnz; A->is_csr = is_csr;
    A->v = malloc(nnz * sizeof(float)); A->idx = malloc(nnz * sizeof(int));
    A->ptr = calloc((out_dim + 1), sizeof(int));
    
    int cur = 0, rem = nnz % out_dim, base = nnz / out_dim;
    for (int i = 0; i < out_dim; i++) {
        int n_elem = base + (i < rem ? 1 : 0);
        for (int p = 0; p < n_elem; p++) {
            A->v[cur] = (rand() % 100) / 10.0f;
            A->idx[cur++] = p * (is_csr ? k : m) / n_elem;
        }
        A->ptr[i + 1] = cur;
    }
}

int main(int argc, char *argv[]) {
    if (argc != 4) {
        printf("Usage: %s <m_rows_A> <k_cols_A_rows_B> <n_cols_B>\n", argv[0]);
        return 1;
    }
    
    int m = atoi(argv[1]), k_dim = atoi(argv[2]), k_cols = atoi(argv[3]);
    
    // Memory allocation
    float *B = malloc(k_dim * k_cols * sizeof(float));
    float *C_csr_seq = malloc(m * k_cols * sizeof(float));
    float *C_csr_omp = malloc(m * k_cols * sizeof(float));
    float *C_csc_seq = malloc(m * k_cols * sizeof(float));
    float *C_mkl_csr = malloc(m * k_cols * sizeof(float));
    float *C_mkl_csc = malloc(m * k_cols * sizeof(float));
    
    for (int i = 0; i < k_dim * k_cols; i++) B[i] = (float)(rand() % 10) / 10.0f;

    SpMat A_csr, A_csc;
    gen_sp(&A_csr, m, k_dim, 1);
    gen_sp(&A_csc, m, k_dim, 0);

    // Custom Benchmarks
    double start = omp_get_wtime(); spmm(&A_csr, B, C_csr_seq, k_cols, 0); 
    printf("CSR (Seq): %f s\n", omp_get_wtime() - start);
    
    start = omp_get_wtime(); spmm(&A_csr, B, C_csr_omp, k_cols, 1); 
    printf("CSR (OMP): %f s\n", omp_get_wtime() - start);

    start = omp_get_wtime(); spmm(&A_csc, B, C_csc_seq, k_cols, 0); 
    printf("CSC (Seq): %f s\n", omp_get_wtime() - start);

    // Intel oneMKL Benchmarks
    sparse_matrix_t mkl_csr, mkl_csc;
    struct matrix_descr descr = {SPARSE_MATRIX_TYPE_GENERAL, 0, 0};
    
    mkl_sparse_s_create_csr(&mkl_csr, SPARSE_INDEX_BASE_ZERO, m, k_dim, 
                            A_csr.ptr, A_csr.ptr+1, A_csr.idx, A_csr.v);
    start = omp_get_wtime();
    mkl_sparse_s_mm(SPARSE_OPERATION_NON_TRANSPOSE, 1.0f, mkl_csr, descr, 
                    SPARSE_LAYOUT_ROW_MAJOR, B, k_cols, k_cols, 0.0f, C_mkl_csr, k_cols);
    printf("oneMKL CSR: %f s\n", omp_get_wtime() - start);

    mkl_sparse_s_create_csc(&mkl_csc, SPARSE_INDEX_BASE_ZERO, m, k_dim, 
                            A_csc.ptr, A_csc.ptr+1, A_csc.idx, A_csc.v);
    start = omp_get_wtime();
    mkl_sparse_s_mm(SPARSE_OPERATION_NON_TRANSPOSE, 1.0f, mkl_csc, descr, 
                    SPARSE_LAYOUT_ROW_MAJOR, B, k_cols, k_cols, 0.0f, C_mkl_csc, k_cols);
    printf("oneMKL CSC: %f s\n", omp_get_wtime() - start);

    // Correctness Verification
    int seq_pass = 1, omp_pass = 1, csc_pass = 1;
    for (int i = 0; i < m * k_cols; i++) {
        if (ABS(C_csr_seq[i] - C_mkl_csr[i]) > eps) seq_pass = 0;
        if (ABS(C_csr_omp[i] - C_mkl_csr[i]) > eps) omp_pass = 0;
        if (ABS(C_csc_seq[i] - C_mkl_csc[i]) > eps) csc_pass = 0;
    }

    printf("\n--- Verification ---\n");
    printf("CSR Sequential: %s\n", seq_pass ? "PASSED" : "FAILED");
    printf("CSR OpenMP:     %s\n", omp_pass ? "PASSED" : "FAILED");
    printf("CSC Sequential: %s\n", csc_pass ? "PASSED" : "FAILED");

    // Cleanup
    mkl_sparse_destroy(mkl_csr);
    mkl_sparse_destroy(mkl_csc);
    free(A_csr.v); free(A_csr.idx); free(A_csr.ptr);
    free(A_csc.v); free(A_csc.idx); free(A_csc.ptr);
    free(B); free(C_csr_seq); free(C_csr_omp); free(C_csc_seq);
    free(C_mkl_csr); free(C_mkl_csc);

    return 0;
}
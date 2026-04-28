#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <mkl.h>

#define ABS(x) ((x) < 0.0f ? -(x) : (x))
#define eps 1e-2

// ---------------------------------------------------------
// 1. Data Structures & "Classes"
// ---------------------------------------------------------
typedef struct DenseMatrix DenseMatrix;
typedef struct SparseMatrix SparseMatrix;

typedef enum {
    FORMAT_CSR_SEQ,
    FORMAT_CSR_OMP,
    FORMAT_CSC
} SparseFormat;

struct DenseMatrix {
    int rows;
    int cols;
    float *values;
    
    float (*get_val)(const DenseMatrix *self, int r, int c);
    int   (*get_idx)(const DenseMatrix *self, int r, int c);
};

static int dense_get_idx_impl(const DenseMatrix *self, int r, int c) {
    return self->cols * r + c; 
}

static float dense_get_val_impl(const DenseMatrix *self, int r, int c) {
    return self->values[self->get_idx(self, r, c)];
}

void init_dense_matrix(DenseMatrix *mat, int rows, int cols, float *values) {
    mat->rows = rows;
    mat->cols = cols;
    mat->values = values;
    mat->get_val = dense_get_val_impl;
    mat->get_idx = dense_get_idx_impl;
}

struct SparseMatrix {
    int rows;
    int cols;
    int nnz;
    SparseFormat fmt;
    float *values;
    int *indices;
    int *indptr;
    
    void (*spmm)(const SparseMatrix *self, const DenseMatrix *B, DenseMatrix *C);
};

// ---------------------------------------------------------
// 2. SpMM Implementations
// ---------------------------------------------------------

static void spmm_csr_seq_impl(const SparseMatrix *self, const DenseMatrix *B, DenseMatrix *C) {
    for (int i = 0; i < self->rows; i++) {
        for (int j = 0; j < B->cols; j++) C->values[C->get_idx(C, i, j)] = 0.0f;
        
        for (int p = self->indptr[i]; p < self->indptr[i+1]; p++) {
            int k = self->indices[p];
            float val = self->values[p];
            for (int j = 0; j < B->cols; j++) {
                C->values[C->get_idx(C, i, j)] += val * B->get_val(B, k, j);
            }
        }
    }
}

static void spmm_csr_omp_impl(const SparseMatrix *self, const DenseMatrix *B, DenseMatrix *C) {
    #pragma omp parallel for
    for (int i = 0; i < self->rows; i++) {
        for (int j = 0; j < B->cols; j++) C->values[C->get_idx(C, i, j)] = 0.0f;
        
        for (int p = self->indptr[i]; p < self->indptr[i+1]; p++) {
            int k = self->indices[p];
            float val = self->values[p];
            for (int j = 0; j < B->cols; j++) {
                C->values[C->get_idx(C, i, j)] += val * B->get_val(B, k, j);
            }
        }
    }
}

static void spmm_csc_impl(const SparseMatrix *self, const DenseMatrix *B, DenseMatrix *C) {
    for (int i = 0; i < C->rows * C->cols; i++) C->values[i] = 0.0f;

    for (int k = 0; k < self->cols; k++) {
        for (int p = self->indptr[k]; p < self->indptr[k+1]; p++) {
            int i = self->indices[p];
            float val = self->values[p];
            for (int j = 0; j < B->cols; j++) {
                C->values[C->get_idx(C, i, j)] += val * B->get_val(B, k, j);
            }
        }
    }
}

void init_sparse_matrix(SparseMatrix *mat, int rows, int cols, int nnz, 
                        SparseFormat fmt, float *values, int *indices, int *indptr) {
    mat->rows = rows;
    mat->cols = cols;
    mat->nnz = nnz;
    mat->fmt = fmt;
    mat->values = values;
    mat->indices = indices;
    mat->indptr = indptr;
    
    if (fmt == FORMAT_CSR_SEQ) mat->spmm = spmm_csr_seq_impl;
    else if (fmt == FORMAT_CSR_OMP) mat->spmm = spmm_csr_omp_impl;
    else mat->spmm = spmm_csc_impl;
}

// ---------------------------------------------------------
// 3. Helper: Generate Random Matrices (NNZ = 32 * m)
// ---------------------------------------------------------
void generate_random_csr_arrays(int m, int k, int total_nnz, 
                                float **val, int **col_idx, int **row_ptr) {
    *val = (float*)malloc(total_nnz * sizeof(float));
    *col_idx = (int*)malloc(total_nnz * sizeof(int));
    *row_ptr = (int*)malloc((m + 1) * sizeof(int));
    
    int current_nnz = 0;
    (*row_ptr)[0] = 0;
    
    int nnz_per_row = total_nnz / m; 
    
    for (int i = 0; i < m; i++) {
        for (int p = 0; p < nnz_per_row; p++) {
            (*val)[current_nnz] = (float)(rand() % 100) / 10.0f;
            (*col_idx)[current_nnz] = (p * k) / nnz_per_row; 
            current_nnz++;
        }
        (*row_ptr)[i + 1] = current_nnz;
    }
}

void generate_random_csc_arrays(int m, int k, int total_nnz, 
                                float **val, int **row_idx, int **col_ptr) {
    *val = (float*)malloc(total_nnz * sizeof(float));
    *row_idx = (int*)malloc(total_nnz * sizeof(int));
    *col_ptr = (int*)malloc((k + 1) * sizeof(int));
    
    int current_nnz = 0;
    (*col_ptr)[0] = 0;
    
    int base_nnz_per_col = total_nnz / k;
    int remainder = total_nnz % k;
    
    for (int j = 0; j < k; j++) {
        int nnz_for_this_col = base_nnz_per_col + (j < remainder ? 1 : 0);
        for (int p = 0; p < nnz_for_this_col; p++) {
            (*val)[current_nnz] = (float)(rand() % 100) / 10.0f;
            (*row_idx)[current_nnz] = (p * m) / nnz_for_this_col; 
            current_nnz++;
        }
        (*col_ptr)[j + 1] = current_nnz;
    }
}

// ---------------------------------------------------------
// 4. Main Benchmark Routine
// ---------------------------------------------------------
int main(int argc, char *argv[]) {
    if (argc != 4) {
        printf("Usage: %s <m_rows_A> <k_cols_A_rows_B> <n_cols_B>\n", argv[0]);
        return 1;
    }

    int m = atoi(argv[1]);
    int k_dim = atoi(argv[2]);
    int k_cols = atoi(argv[3]);
    int total_nnz = 32 * m; 

    srand(42);
    printf("Benchmarking A(%d x %d) * B(%d x %d) -> C(%d x %d)\n", m, k_dim, k_dim, k_cols, m, k_cols);
    printf("Total NNZ: %d\n\n", total_nnz);

    // 1. Setup Dense B
    float *b_vals = (float*)malloc(k_dim * k_cols * sizeof(float));
    for (int i = 0; i < k_dim * k_cols; i++) b_vals[i] = (float)(rand() % 10) / 10.0f;
    DenseMatrix B;
    init_dense_matrix(&B, k_dim, k_cols, b_vals);

    // ---------------------------------------------------------
    // CSR Benchmarking Block
    // ---------------------------------------------------------
    printf("--- CSR Format Testing ---\n");
    float *csr_val; int *csr_col_idx, *csr_row_ptr;
    generate_random_csr_arrays(m, k_dim, total_nnz, &csr_val, &csr_col_idx, &csr_row_ptr);

    float *c_seq_vals = (float*)malloc(m * k_cols * sizeof(float));
    float *c_omp_vals = (float*)malloc(m * k_cols * sizeof(float));
    float *c_mkl_csr_vals = (float*)malloc(m * k_cols * sizeof(float));

    DenseMatrix C_seq, C_omp;
    init_dense_matrix(&C_seq, m, k_cols, c_seq_vals);
    init_dense_matrix(&C_omp, m, k_cols, c_omp_vals);
    SparseMatrix A_seq, A_omp;
    init_sparse_matrix(&A_seq, m, k_dim, total_nnz, FORMAT_CSR_SEQ, csr_val, csr_col_idx, csr_row_ptr);
    init_sparse_matrix(&A_omp, m, k_dim, total_nnz, FORMAT_CSR_OMP, csr_val, csr_col_idx, csr_row_ptr);

    double start = omp_get_wtime(); A_seq.spmm(&A_seq, &B, &C_seq); double t_seq = omp_get_wtime() - start;
    start = omp_get_wtime(); A_omp.spmm(&A_omp, &B, &C_omp); double t_omp = omp_get_wtime() - start;

    sparse_matrix_t mkl_A_csr;
    mkl_sparse_s_create_csr(&mkl_A_csr, SPARSE_INDEX_BASE_ZERO, m, k_dim, csr_row_ptr, csr_row_ptr+1, csr_col_idx, csr_val);
    struct matrix_descr descr = {SPARSE_MATRIX_TYPE_GENERAL, 0, 0};
    
    mkl_sparse_s_mm(SPARSE_OPERATION_NON_TRANSPOSE, 1.0f, mkl_A_csr, descr, SPARSE_LAYOUT_ROW_MAJOR, B.values, k_cols, k_cols, 0.0f, c_mkl_csr_vals, k_cols);
    start = omp_get_wtime();
    mkl_sparse_s_mm(SPARSE_OPERATION_NON_TRANSPOSE, 1.0f, mkl_A_csr, descr, SPARSE_LAYOUT_ROW_MAJOR, B.values, k_cols, k_cols, 0.0f, c_mkl_csr_vals, k_cols);
    double t_mkl_csr = omp_get_wtime() - start;

    printf("Custom CSR (Seq): %f s\n", t_seq);
    printf("Custom CSR (OMP): %f s\n", t_omp);
    printf("Intel oneMKL CSR: %f s\n", t_mkl_csr);

    // ---------------------------------------------------------
    // CSC Benchmarking Block
    // ---------------------------------------------------------
    printf("\n--- CSC Format Testing ---\n");
    float *csc_val; int *csc_row_idx, *csc_col_ptr;
    generate_random_csc_arrays(m, k_dim, total_nnz, &csc_val, &csc_row_idx, &csc_col_ptr);

    float *c_csc_vals = (float*)malloc(m * k_cols * sizeof(float));
    float *c_mkl_csc_vals = (float*)malloc(m * k_cols * sizeof(float));

    DenseMatrix C_csc;
    init_dense_matrix(&C_csc, m, k_cols, c_csc_vals);
    SparseMatrix A_csc;
    init_sparse_matrix(&A_csc, m, k_dim, total_nnz, FORMAT_CSC, csc_val, csc_row_idx, csc_col_ptr);

    start = omp_get_wtime(); A_csc.spmm(&A_csc, &B, &C_csc); double t_csc = omp_get_wtime() - start;

    sparse_matrix_t mkl_A_csc;
    mkl_sparse_s_create_csc(&mkl_A_csc, SPARSE_INDEX_BASE_ZERO, m, k_dim, csc_col_ptr, csc_col_ptr+1, csc_row_idx, csc_val);
    
    mkl_sparse_s_mm(SPARSE_OPERATION_NON_TRANSPOSE, 1.0f, mkl_A_csc, descr, SPARSE_LAYOUT_ROW_MAJOR, B.values, k_cols, k_cols, 0.0f, c_mkl_csc_vals, k_cols);
    start = omp_get_wtime();
    mkl_sparse_s_mm(SPARSE_OPERATION_NON_TRANSPOSE, 1.0f, mkl_A_csc, descr, SPARSE_LAYOUT_ROW_MAJOR, B.values, k_cols, k_cols, 0.0f, c_mkl_csc_vals, k_cols);
    double t_mkl_csc = omp_get_wtime() - start;

    printf("Custom CSC (Seq): %f s\n", t_csc);
    printf("Intel oneMKL CSC: %f s\n", t_mkl_csc);

    // ---------------------------------------------------------
    // Verification
    // ---------------------------------------------------------
    int seq_pass = 1, omp_pass = 1, csc_pass = 1;
    for (int i = 0; i < m * k_cols; i++) {
        if (ABS(C_seq.values[i] - c_mkl_csr_vals[i]) > eps) seq_pass = 0;
        if (ABS(C_omp.values[i] - c_mkl_csr_vals[i]) > eps) omp_pass = 0;
        if (ABS(C_csc.values[i] - c_mkl_csc_vals[i]) > eps) csc_pass = 0;
    }

    printf("\n--- Verification ---\n");
    printf("CSR Sequential: %s\n", seq_pass ? "PASSED" : "FAILED");
    printf("CSR OpenMP:     %s\n", omp_pass ? "PASSED" : "FAILED");
    printf("CSC Sequential: %s\n", csc_pass ? "PASSED" : "FAILED");

    // Cleanup
    mkl_sparse_destroy(mkl_A_csr);
    mkl_sparse_destroy(mkl_A_csc);
    free(csr_val); free(csr_col_idx); free(csr_row_ptr);
    free(csc_val); free(csc_row_idx); free(csc_col_ptr);
    free(b_vals); free(c_seq_vals); free(c_omp_vals); free(c_csc_vals);
    free(c_mkl_csr_vals); free(c_mkl_csc_vals);

    return 0;
}
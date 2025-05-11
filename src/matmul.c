#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <omp.h>
#include "matmul.h"
#include "utils.h"

// Forward declaration of the recursive helper for parallel D&C
void _dac_task(long *A, long *B, long *C, int n, int base_thresh);

int main(int argc, char **argv) {
    // Command‐line argument variables
    char *input = NULL, *method = NULL;
    int base_thresh = 64;    // default base case for recursion
    int threads = 0;         // default to using OMP_NUM_THREADS
    int opt;

    // Parse flags: -i input, -m method, -b base threshold, -t thread count
    while ((opt = getopt(argc, argv, "i:m:b:t:")) != -1) {
        switch (opt) {
            case 'i': input       = optarg;         break;
            case 'm': method      = optarg;         break;
            case 'b': base_thresh = atoi(optarg);   break;
            case 't': threads     = atoi(optarg);   break;
            default:
                fprintf(stderr,
                    "Usage: %s -i <input> -m <method> [-b <base>] [-t <threads>]\n",
                    argv[0]
                );
                exit(EXIT_FAILURE);
        }
    }

    // Validate required flags
    if (!input || !method) {
        fprintf(stderr, "Error: Must specify both -i and -m flags\n");
        exit(EXIT_FAILURE);
    }

    // If threads flag set, override OpenMP thread count
    if (threads > 0) {
        omp_set_num_threads(threads);
    }

    // Read matrices A and B from file; get dimension n
    long *A, *B, *C;
    int n;
    read_input(input, &A, &B, &n);

    // Allocate result matrix C and zero‐initialize
    C = malloc(sizeof(long) * n * n);
    memset(C, 0, sizeof(long) * n * n);

    // Start timing (exclude I/O)
    double t0 = omp_get_wtime();

    // Dispatch to the selected multiplication method
    if      (strcmp(method, "Sequential") == 0) {
        mult_seq(A, B, C, n);
    }
    else if (strcmp(method, "SequentialP") == 0) {
        mult_seq_par(A, B, C, n);
    }
    else if (strcmp(method, "StraightDivAndConq") == 0) {
        mult_dac(A, B, C, n, base_thresh);
    }
    else if (strcmp(method, "StraightDivAndConqP") == 0) {
        mult_dac_par(A, B, C, n, base_thresh);
    }
    else if (strcmp(method, "StrassenDivAndConq") == 0) {
        mult_strassen(A, B, C, n, base_thresh);
    }
    else {
        fprintf(stderr, "Error: Unknown method '%s'\n", method);
        exit(EXIT_FAILURE);
    }

    // Stop timing and calculate elapsed seconds
    double elapsed = omp_get_wtime() - t0;

    // Build output filenames based on input name, dimension, and method
    char outM[256], outI[256];

    // Extract filename without path or “.txt”
    char *base = strrchr(input, '/');
    base = base ? base + 1 : input;
    char name[256];
    strncpy(name, base, strlen(base) - 4);
    name[strlen(base) - 4] = '\0';

    snprintf(outM, sizeof(outM), "outputs/%s_%d_output_%s.txt", name, n, method);
    snprintf(outI, sizeof(outI), "outputs/%s_%d_info_%s.txt", name, n, method);

    // Write resulting matrix and timing info to files
    write_matrix(outM, C, n);
    write_time(outI, elapsed);

    // Clean up
    free(A);
    free(B);
    free(C);
    return EXIT_SUCCESS;
}

//---------------- Sequential Methods ----------------//

/**
 * Basic triple‐loop matrix multiplication (i,k,j order).
 * C[i][j] += A[i][k] * B[k][j]
 */
void mult_seq(long *A, long *B, long *C, int n) {
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) {
            long aik = A[i*n + k];
            for (int j = 0; j < n; j++) {
                C[i*n + j] += aik * B[k*n + j];
            }
        }
    }
}

/**
 * Parallelized triple‐loop using OpenMP
 * Each (i,j) pair is computed in parallel threads
 */
void mult_seq_par(long *A, long *B, long *C, int n) {
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) {
            long aik = A[i*n + k];
            for (int j = 0; j < n; j++) {
                C[i*n + j] += aik * B[k*n + j];
            }
        }
    }
}

//--------------- Divide & Conquer Methods ---------------//

/**
 * Sequential D&C: split each matrix into four n/2 × n/2 subblocks,
 * recursively multiply and accumulate.
 */
void mult_dac(long *A, long *B, long *C, int n, int base_thresh) {
    // Base case: fall back to simple multiplication
    if (n <= base_thresh) {
        mult_seq(A, B, C, n);
        return;
    }

    int m = n / 2;
    size_t sz = (size_t)m * m;

    // Compute pointers to each subblock
    long *A11 = A;
    long *A12 = A + m;
    long *A21 = A + m*n;
    long *A22 = A + m*n + m;
    long *B11 = B;
    long *B12 = B + m;
    long *B21 = B + m*n;
    long *B22 = B + m*n + m;
    long *C11 = C;
    long *C12 = C + m;
    long *C21 = C + m*n;
    long *C22 = C + m*n + m;

    // Temporary buffer for partial results
    long *T = calloc(sz, sizeof(long));

    // C11 = A11*B11 + A12*B21
    mult_dac(A11, B11, C11, m, base_thresh);
    mult_dac(A12, B21, T,   m, base_thresh);
    for (size_t i = 0; i < sz; i++) C11[i] += T[i];

    // C12 = A11*B12 + A12*B22
    mult_dac(A11, B12, C12, m, base_thresh);
    mult_dac(A12, B22, T,   m, base_thresh);
    for (size_t i = 0; i < sz; i++) C12[i] += T[i];

    // C21 = A21*B11 + A22*B21
    mult_dac(A21, B11, C21, m, base_thresh);
    mult_dac(A22, B21, T,   m, base_thresh);
    for (size_t i = 0; i < sz; i++) C21[i] += T[i];

    // C22 = A21*B12 + A22*B22
    mult_dac(A21, B12, C22, m, base_thresh);
    mult_dac(A22, B22, T,   m, base_thresh);
    for (size_t i = 0; i < sz; i++) C22[i] += T[i];

    free(T);
}

/**
 * Helper for parallel D&C: spawns tasks recursively.
 */
void _dac_task(long *A, long *B, long *C, int n, int base_thresh) {
    if (n <= base_thresh) {
        mult_seq(A, B, C, n);
        return;
    }
    int m = n / 2;
    size_t sz = (size_t)m * m;

    // Pointers as before
    long *A11 = A, *A12 = A + m, *A21 = A + m*n, *A22 = A + m*n + m;
    long *B11 = B, *B12 = B + m, *B21 = B + m*n, *B22 = B + m*n + m;
    long *C11 = C, *C12 = C + m, *C21 = C + m*n, *C22 = C + m*n + m;

    // Two temporaries to hold partial results
    long *T1 = calloc(sz, sizeof(long));
    long *T2 = calloc(sz, sizeof(long));

    // Launch 8 recursive tasks
    #pragma omp task
    _dac_task(A11, B11, C11, m, base_thresh);
    #pragma omp task
    _dac_task(A12, B21, T1,  m, base_thresh);
    #pragma omp task
    _dac_task(A11, B12, C12, m, base_thresh);
    #pragma omp task
    _dac_task(A12, B22, T2,  m, base_thresh);
    #pragma omp task
    _dac_task(A21, B11, C21, m, base_thresh);
    #pragma omp task
    _dac_task(A22, B21, T1,  m, base_thresh);  // reuse T1
    #pragma omp task
    _dac_task(A21, B12, C22, m, base_thresh);
    #pragma omp task
    _dac_task(A22, B22, T2,  m, base_thresh);  // reuse T2
    #pragma omp taskwait

    // Accumulate partials
    for (size_t i = 0; i < sz; i++) C11[i] += T1[i];
    for (size_t i = 0; i < sz; i++) C12[i] += T2[i];
    for (size_t i = 0; i < sz; i++) C21[i] += T1[i];
    for (size_t i = 0; i < sz; i++) C22[i] += T2[i];

    free(T1);
    free(T2);
}

/**
 * Parallel D&C entry: starts the OpenMP parallel region
 */
void mult_dac_par(long *A, long *B, long *C, int n, int base_thresh) {
    #pragma omp parallel
    #pragma omp single
    _dac_task(A, B, C, n, base_thresh);
}

//--------------------- Strassen’s ----------------------//

/**
 * Strassen’s algorithm: 7 recursive multiplies + combines.
 */
void mult_strassen(long *A, long *B, long *C, int n, int base_thresh) {
    if (n <= base_thresh) {
        mult_seq(A, B, C, n);
        return;
    }
    int m = n / 2;
    size_t sz = (size_t)m * m;

    // Allocate temporary arrays S1,S2 for sums/diffs, and P for product
    long *S1 = malloc(sz * sizeof(long));
    long *S2 = malloc(sz * sizeof(long));
    long *P  = malloc(sz * sizeof(long));

    // Define subblocks as before
    long *A11=A,      *A12=A + m;
    long *A21=A + m*n, *A22=A + m*n + m;
    long *B11=B,      *B12=B + m;
    long *B21=B + m*n, *B22=B + m*n + m;
    long *C11=C,      *C12=C + m;
    long *C21=C + m*n, *C22=C + m*n + m;

    // 1) P1 = (A11+A22)*(B11+B22)
    for (size_t i = 0; i < sz; i++) S1[i] = A11[i] + A22[i];
    for (size_t i = 0; i < sz; i++) S2[i] = B11[i] + B22[i];
    mult_strassen(S1, S2, P, m, base_thresh);
    for (size_t i = 0; i < sz; i++) C11[i] = P[i];
    for (size_t i = 0; i < sz; i++) C22[i] = P[i];

    // 2) P2 = (A21 + A22)*B11
    for(size_t i=0;i<sz;i++) S1[i] = A21[i] + A22[i];
    mult_strassen(S1,B11,P,m,base_thresh);
    for(size_t i=0;i<sz;i++) C21[i] = P[i];
    for(size_t i=0;i<sz;i++) C11[i] += P[i];

    // 3) P3 = A11*(B12 - B22)
    for(size_t i=0;i<sz;i++) S2[i] = B12[i] - B22[i];
    mult_strassen(A11,S2,P,m,base_thresh);
    for(size_t i=0;i<sz;i++) C12[i] = P[i];
    for(size_t i=0;i<sz;i++) C11[i] += P[i];

    // 4) P4 = A22*(B21 - B11)
    for(size_t i=0;i<sz;i++) S2[i] = B21[i] - B11[i];
    mult_strassen(A22,S2,P,m,base_thresh);
    for(size_t i=0;i<sz;i++) C22[i] += P[i];
    for(size_t i=0;i<sz;i++) C12[i] -= P[i];

    // 5) P5 = (A11 + A12)*B22
    for(size_t i=0;i<sz;i++) S1[i] = A11[i] + A12[i];
    mult_strassen(S1,B22,P,m,base_thresh);
    for(size_t i=0;i<sz;i++) C12[i] += P[i];
    for(size_t i=0;i<sz;i++) C22[i] += P[i];

    // 6) P6 = (A21 - A11)*(B11 + B12)
    for(size_t i=0;i<sz;i++) S1[i] = A21[i] - A11[i];
    for(size_t i=0;i<sz;i++) S2[i] = B11[i] + B12[i];
    mult_strassen(S1,S2,P,m,base_thresh);
    for(size_t i=0;i<sz;i++) C11[i] -= P[i];
    for(size_t i=0;i<sz;i++) C21[i] += P[i];

    // 7) P7 = (A12 - A22)*(B21 + B22)
    for(size_t i=0;i<sz;i++) S1[i] = A12[i] - A22[i];
    for(size_t i=0;i<sz;i++) S2[i] = B21[i] + B22[i];
    mult_strassen(S1,S2,P,m,base_thresh);
    for(size_t i=0;i<sz;i++) C12[i] -= P[i];
    for(size_t i=0;i<sz;i++) C22[i] += P[i];

    free(S1);
    free(S2);
    free(P);
}
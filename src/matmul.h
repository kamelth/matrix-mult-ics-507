#ifndef MATMUL_H
#define MATMUL_H

// Sequential & parallel straightforward
void mult_seq(long *A, long *B, long *C, int n);
void mult_seq_par(long *A, long *B, long *C, int n);

// Divide & conquer (threshold = base case size)
void mult_dac(long *A, long *B, long *C, int n, int base_thresh);
void mult_dac_par(long *A, long *B, long *C, int n, int base_thresh);

// Strassen’s algorithm (threshold = when to switch to simple multiply)
void mult_strassen(long *A, long *B, long *C, int n, int base_thresh);

#endif

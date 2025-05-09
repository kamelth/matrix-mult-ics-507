#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>

// Read two n×n matrices from 'fname':
//  1st line = n
//  next n lines = A (row-major), then n lines = B
void read_input(const char *fname, long **A, long **B, int *n);

// Write an n×n matrix M to 'fname' (space-separated rows)
void write_matrix(const char *fname, long *M, int n);

// Write elapsed time (secs) as hh:mm:ss into 'fname'
void write_time(const char *fname, double secs);

// Utility: format seconds → "HH:MM:SS"
void format_hhmmss(double secs, char *buf);

#endif

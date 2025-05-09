// utils.c
#include "utils.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdio.h>

void format_hhmmss(double secs, char *buf) {
    int h = (int)secs / 3600;
    int m = ((int)secs % 3600) / 60;
    int s = (int)secs % 60;
    sprintf(buf, "%02d:%02d:%02d", h, m, s);
}

void read_input(const char *fname, long **A, long **B, int *n) {
    FILE *f = fopen(fname, "r");
    if (!f) { perror("fopen"); exit(1); }
    fscanf(f, "%d", n);
    int N = *n;
    *A = malloc(sizeof(long)*N*N);
    *B = malloc(sizeof(long)*N*N);
    for (int i = 0; i < N*N; i++) {
        fscanf(f, "%ld", &(*A)[i]);
    }
    for (int i = 0; i < N*N; i++) {
        fscanf(f, "%ld", &(*B)[i]);
    }
    fclose(f);
}

void write_matrix(const char *fname, long *M, int n) {
    FILE *f = fopen(fname, "w");
    if (!f) { perror("fopen"); exit(1); }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fprintf(f, "%ld%s", M[i*n + j], (j+1==n ? "\n" : " "));
        }
    }
    fclose(f);
}

void write_time(const char *fname, double secs) {
    char buf[16];
    format_hhmmss(secs, buf);
    FILE *f = fopen(fname, "w");
    if (!f) { perror("fopen"); exit(1); }
    fprintf(f, "%s\n", buf);
    fclose(f);
}

#include <stdio.h>

void
matmult_nat(int m, int n, int k, double **A, double **B, double **C) {
    
    int i, j, l;

    for(i = 0; i < m; i++){
        for(j = 0; j < n; j++){
            double tmp = 0.0;
            for(l = 0; l < k; l++){
                tmp = tmp + A[i][l] * B[l][j];
            }
            C[i][j] = tmp;
        }
    }
}

#include <stdio.h>
#include <sunperf.h>

void
matmult_nat(int m, int n, int k, double **A, double **B, double **C) {
    
    int i, j, l;
    

    for(j = 0; j < n; j++){
        for(i = 0; i < m; i++){
	    double tmp = 0.0;
            for(l = 0; l < k; l++){
                tmp += A[i][l] * B[l][j];
            }
            C[i][j] = tmp;
        }
    }
}


void
matmult_lib(int m, int n, int k, double **A, double **B, double **C) {
	dgemm('N', 'N', m ,m , n, 1.0, *A, k, *B, n, 0.0, *C, n);
}


// MNK

void
matmult_mnk(int m, int n, int k, double **A, double **B, double **C) {
    
    int i, j, l;

    for(i = 0; i < m; i++){
        for(j = 0; j < n; j++){
		C[i][j] = 0;
    	}
    }


    for(i = 0; i < m; i++){
        for(j = 0; j < n; j++){
            for(l = 0; l < k; l++){
                C[i][j] += A[i][l] * B[l][j];
            }
        }
    }
}


// NMK

void
matmult_nmk(int m, int n, int k, double **A, double **B, double **C) {
    
    int i, j, l;

    for(i = 0; i < m; i++){
        for(j = 0; j < n; j++){
		C[i][j] = 0;
    	}
    }

    for(j = 0; j < n; j++){
        for(i = 0; i < m; i++){
            for(l = 0; l < k; l++){
                C[i][j] += A[i][l] * B[l][j];
            }
        }
    }
}


// KMN

void
matmult_kmn(int m, int n, int k, double **A, double **B, double **C) {
    
    int i, j, l;

    for(i = 0; i < m; i++){
        for(j = 0; j < n; j++){
		C[i][j] = 0;
    	}
    }

    for(l = 0; l < k; l++){
        for(i = 0; i < m; i++){
            for(j = 0; j < n; j++){
                C[i][j] += A[i][l] * B[l][j];
            }
        }
    }
}


// MKN

void
matmult_mkn(int m, int n, int k, double **A, double **B, double **C) {
    
    int i, j, l;

    for(i = 0; i < m; i++){
        for(j = 0; j < n; j++){
		C[i][j] = 0;
    	}
    }


    for(i = 0; i < m; i++){
        for(l = 0; l < k; l++){
            for(j = 0; j < n; j++){
                C[i][j] += A[i][l] * B[l][j];
            }
        }
    }
}


// NKM

void
matmult_nkm(int m, int n, int k, double **A, double **B, double **C) {
    
    int i, j, l;

    for(i = 0; i < m; i++){
        for(j = 0; j < n; j++){
		C[i][j] = 0;
    	}
    }

    for(j = 0; j < n; j++){
        for(l = 0; l < k; l++){
            for(i = 0; i < m; i++){
                C[i][j] += A[i][l] * B[l][j];
            }
        }
    }
}


// KNM

void
matmult_knm(int m, int n, int k, double **A, double **B, double **C) {
    
    int i, j, l;

    for(i = 0; i < m; i++){
        for(j = 0; j < n; j++){
		C[i][j] = 0;
    	}
    }

    for(l = 0; l < k; l++){
        for(j = 0; j < n; j++){
            for(i = 0; i < m; i++){
                C[i][j] += A[i][l] * B[l][j];
            }
        }
    }
}

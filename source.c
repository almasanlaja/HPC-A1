#include <stdio.h>
#include <sunperf.h>
#include <math.h>

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


// Blocks

void 
matmult_blk(int m, int n, int k, double **A, double **B, double **C, int bs) {

    int i, j, l;

    for(i = 0; i < m; i++){
        for(j = 0; j < n; j++){
		C[i][j] = 0;
    	}
    }

    int ns, ms, ks, rns, rms, rks;
	
    ns = n/bs;
    ms = m/bs;
    ks = k/bs;

    rns = n%bs;
    rms = m%bs;
    rks = k%bs;

    int is, js, ls;

    for (is=0; is < ms; is++){
	for(js=0; js < ns; js++){
		for(ls=0; ls < ks; ls++){
			for(i = is*bs; i < (is+1)*bs; i++){
				for(j = js*bs; j < (js+1)*bs; j++){
				    double tmp1 = 0.0;
				    for(l = ls*bs; l < (ls+1)*bs; l++){
					tmp1 += A[i][l] * B[l][j];
				    }
				    C[i][j] += tmp1;
				}
    			}
		}
		
		for(i = is*bs; i < (is+1)*bs; i++){
			for(j = js*bs; j < (js+1)*bs; j++){
			    double tmp2 = 0.0;
			    for(l = k-rks; l < k; l++){
				tmp2 += A[i][l] * B[l][j];
			    }
			    C[i][j] += tmp2;
			}
		}
	}
	for(ls=0; ls < ks; ls++){
		for(i = is*bs; i < (is+1)*bs; i++){
			for(j = n-rns; j < n; j++){
			    double tmp3 = 0.0;
			    for(l = ls*bs; l < (ls+1)*bs; l++){
				tmp3 += A[i][l] * B[l][j];
			    }
			    C[i][j] += tmp3;
			}
		}
	}
	
	for(i = is*bs; i < (is+1)*bs; i++){
		for(j = n-rns; j < n; j++){
		    double tmp4 = 0.0;
		    for(l = k-rks; l < k; l++){
			tmp4 += A[i][l] * B[l][j];
		    }
		    C[i][j] += tmp4;
		}
	}

    }


    for(js=0; js < ns; js++){
		for(ls=0; ls < ks; ls++){
			for(i = m-rms; i < m; i++){
				for(j = js*bs; j < (js+1)*bs; j++){
				    double tmp5 = 0.0;
				    for(l = ls*bs; l < (ls+1)*bs; l++){
					tmp5 += A[i][l] * B[l][j];
				    }
				    C[i][j] += tmp5;
				}
			}
		}
		
		for(i = m-rms; i < m; i++){
			for(j = js*bs; j < (js+1)*bs; j++){
			    double tmp6 = 0.0;
			    for(l = k-rks; l < k; l++){
				tmp6 += A[i][l] * B[l][j];
			    }
			    C[i][j] += tmp6;
			}
		}
	}
	for(ls=0; ls < ks; ls++){
		for(i = m-rms; i < m; i++){
			for(j = n-rns; j < n; j++){
			    double tmp7 = 0.0;
			    for(l = ls*bs; l < (ls+1)*bs; l++){
				tmp7 += A[i][l] * B[l][j];
			    }
			    C[i][j] += tmp7;
			}
		}
	}

	for(i = m-rms; i < m; i++){
		for(j = n-rns; j < n; j++){
		    double tmp8 = 0.0;
		    for(l = k-rks; l < k; l++){
			tmp8 += A[i][l] * B[l][j];
		    }
		    C[i][j] += tmp8;
		}
	}
    
    
}

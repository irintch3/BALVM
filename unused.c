/*******************************************************************/


void mu_comp(double *eta, int p, int n, int order, int m, int nknots, double *knots, double *knotsI, double *all_beta, double *all_mu)
{
/* 
m is the degree of freedom of splines and number of components in beta_{kj}
knots and knots is for splines without intercept (but the first coordinates is the intercept function and should be excluded)
number of knots in knotsI is knots-1
order is the order of the spline which is the degree+1
n is the number od observations
*/
 
  int i, j,k, supi;
  double *B, *BI; 
  int *offsets, *offsetsI, *derivs;
  
  B=(double*)calloc(order*n*p, sizeof(double));
  BI=(double*)calloc(order*n*p, sizeof(double));
  
  derivs=(int*)calloc(n*p, sizeof(int));    
  for(i = 0; i < n*p; i++) derivs[i] = 0; 
  
  offsets=(int*)calloc(n*p, sizeof(int));
  offsetsI=(int*)calloc(n*p, sizeof(int));
  
  spline_basis(knotsI, (nknots-1), order, eta, n*p, derivs, n*p, offsetsI, BI);
  print_matrix( "BI: ", order, n, BI, n );
  for (i=0; i<n*p; i++){Rprintf("  offsetsI[%d]=%d\n",i, offsetsI[i]);}
  spline_basis(knots, nknots, order, eta, n*p, derivs, n*p, offsets, B);
  print_matrix( "B: ", order, n, B, n );
  
  
  for(j = 0; j < p; j++) {
     for (k=0; k<=j; k++){
       if (k==0) {
	for (i=0; i<n; i++) {
	  for (supi=0; supi<order; supi++){
	    all_mu[j*(j+1)*n/2+k]=all_mu[j*(j+1)*n/2+k]+BI[j*order*n+i*order+supi]*all_beta[j*(j+1)*m/2+offsetsI[j*n+i]-1+supi];
					  }
			      }
		} else {
		    for (i=0; i<n; i++) {
		      if (offsets[j*n+i]!=0) {
			for (supi=0; supi<order; supi++){
			  all_mu[j*(j+1)*n/2+k*n+i]=all_mu[j*(j+1)*n/2+k*n+i]+B[j*order*n+i*order+supi]*all_beta[ j*(j+1)*m/2+k*m+offsets[j*n+i]-1+supi ];
			}
						  } else {for (supi=0; supi<(order-1); supi++){
						    all_mu[j*(j+1)*n/2+k*n+i]=all_mu[j*(j+1)*n/2+k*n+i]+B[j*order*n+i*order+supi+1]*all_beta[ j*(j+1)*m/2+k*m+offsets[j*n+i]+supi ];
												}
							}
												
					}
			}
			}
			  }
   
  free(B);
  free(BI);
  free(derivs);
  free(offsets);
  free(offsetsI);

  }


SEXP sample_check(SEXP Rx, SEXP Rsize, SEXP Rprob)
{
    int i, n, *perm, size, *num_grid;
    double *prob, *x;
    SEXP ans;

    size=INTEGER(Rsize)[0];
    n=length(Rx);
    Rprintf("n=%d\n", n);
    x = (double *) calloc(n, sizeof(double));
    prob = (double *) calloc(n, sizeof(double));
    perm = (int *) calloc(n, sizeof(int));
    num_grid = (int *) calloc(1, sizeof(int));
    for (i = 0; i < n; i++) {x[i]=REAL(Rx)[i]; prob[i]=REAL(Rprob)[i];}

	    int nc = 0;
	    for (i = 0; i < n; i++) if(n * prob[i] > 0.1) nc++;
	    
	    if (nc > 200)
		walker_ProbSampleReplace(n, prob, perm, 1, num_grid);
	    else
		ProbSampleReplace(n, prob, perm, 1, num_grid);  
Rprintf("num_grid[0]=%d\n", num_grid[0]);	    
 PROTECT(ans = allocVector(REALSXP, 1));
 REAL(ans)[0]=x[num_grid[0]-1];
 UNPROTECT(1);
 free(x);
 free(perm);
 free(num_grid);
 free(prob);
 return ans;
}




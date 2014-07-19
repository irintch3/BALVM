/* 
* here are the functions 
* all_mu_comp 
* comp_etaY_grid
* mu_B_fct
* betaj_post_check
* gibbs_BALVM_naive
* B_mu_check
* removed from fcts_Eta_B.c
* to avoid conflict with function existing in Naive.c
*/



/*
* given eta and all_beta computes all_mu
* m is the degree of freedom of splines and number of components in beta_{kj}
* knots and knots is for splines without intercept (but the first coordinates is the intercept function and 
* should be excluded)
* number of knots in knotsI is knots-1
* order is the order of the spline which is the degree+1
* n is the number od observations
* p is the number variables in initial y and defines the size of eta and all_beta
*/
void all_mu_comp(double *eta, int p, int order, int m, int nknots, double *knots, double *knotsI, 
                double *all_beta, double *all_mu)
{
 
  int i, j,k, supi;
  double *B, *BI, *Ball; 
  int *offsets, *offsetsI, *derivs;
  
  B=(double*)calloc(order*p, sizeof(double));
  BI=(double*)calloc(order*p, sizeof(double));
  Ball=(double*)calloc(m*p, sizeof(double));
  
  derivs=(int*)calloc(p, sizeof(int));    
  for(i = 0; i < p; i++) derivs[i] = 0; 
  
  offsets=(int*)calloc(p, sizeof(int));
  offsetsI=(int*)calloc(p, sizeof(int));
  
  spline_basis(knotsI, (nknots-1), order, eta, p, derivs, p, offsetsI, BI);
  /*print_matrix( "BI: ", order, 1, BI, 1 );*/
  /*for (i=0; i<p; i++){Rprintf("  offsetsI[%d]=%d\n",i, offsetsI[i]);}*/
  spline_basis(knots, nknots, order, eta, p, derivs, p, offsets, B);
  /*print_matrix( "B: ", order*p, 1, B, 1 );
  for (i=0; i<p; i++){Rprintf("  offsets[%d]=%d\n",i, offsets[i]);}*/
  
 for (i=0; i<order; i++){
 /*Rprintf("offsetsI[0]-1+i=%d\n", offsetsI[0]-1+i); */
    Ball[offsetsI[0]+i]=BI[i];
  }
  
   for (j=1; j<p; j++){
    if (offsets[j]!=0) {
      for (i=0; i<order; i++){Ball[m*j+offsets[j]+i-1]=B[j*order+i];  /**/
      /*Rprintf("m*j+offsets[j]+i=%d\n",m*j+offsets[j]+i);*/ }
    } else {for (i=1; i<order; i++){ Ball[m*j+offsets[j]+i-1]=B[j*order+i];/**/
        /*Rprintf("m*j+offsets[j]+i=%d\n",m*j+offsets[j]+i)*/;}}
  		}
 /*print_matrix( "Ball: ", m*p,1, Ball, m*p); */
 
  for (j=0;j<p; j++){
    for (k=0; k<=j; k++){  
      for (supi=k*m; supi<(k+1)*m; supi++) all_mu[j*(j+1)/2+k]+=Ball[supi] * all_beta[j*(j+1)*m/2+supi];
    }
  }
free(B);
free(BI);
free(Ball);
free(offsets);
free(offsetsI);
free(derivs);
}

/*
* This procedure computes and simulates the new eta and ystar
*/
void comp_etaY_grid (int n, int p, int S, int order, int m, int Nknots, 
		     double *knots, double *knotsI,  double *G, double *all_beta, double *all_sigma,
		     double *data, double *eta, double *ystar, double ***Pijs)
{
  double *all_mu, *yi_jp,  *mu_p, *mu_jp, *eta_i, *sigma_jp;
  double /*Pijs[n][p][S], */Pij[n][p], Pjs[p][S], data_mat[n][p], eta_mat[n][p];
  int s, i,j,k, j1;
  
  for (i=0; i<n;i++){for (j=0; j<p; j++) {eta_mat[i][j]=eta[j*n+i]; data_mat[i][j]=data[j*n+i]; Pij[i][j]=0;}}
  for (j=0; j<p;j++){for (s=0; s<S; s++) {Pjs[j][s]=0;}}

  eta_i=(double *) calloc(p, sizeof(double));
  
  
for (i=0; i<n; i++){
  for (j=0; j<p; j++){
    for (s=0; s<S; s++){
	  /*for (k=0; k<p; k++) {eta_i[k]=eta_mat[i][k];}
	  eta_i[j]=G[s];*/
	  all_mu=(double *) calloc(p*(p+1)/2, sizeof(double));
	  mu_p=(double *) calloc(p, sizeof(double));
	  all_mu_comp(eta_i, p, order, m, Nknots, knots, knotsI, all_beta, all_mu);/**/
	  /*print_matrix( "all_mu: ", 1, p*(p+1)/2, all_mu, 1 ); */
	  /*for (j1=0; j1<p; j1++){for (k=0;k<j1+1; k++){ mu_p[j1]=mu_p[j1]+all_mu[j1*(j1+1)/2+k]; }}*/
	  /*print_matrix( "mu_p: ", 1, p, mu_p, 1 ); */
	  /*sigma_jp=(double *) calloc(p-j, sizeof(double));
	  yi_jp=(double *) calloc(p-j, sizeof(double));  
	  mu_jp=(double *) calloc(p-j, sizeof(double));*/
	  /*for (k=0; k<p-j; k++) { yi_jp[k]=data_mat[i][j+k]; sigma_jp[k]=all_sigma[j+k]; mu_jp[k]=mu_p[j+k];}
	  Pijs[i][j][s]=dmvnorm(p-j, sigma_jp, yi_jp, mu_jp);*/
	  /*Rprintf( "i=%d  j=%d  s=%d   ", i+1, j+1, s+1 ); */
	  /*print_matrix( "yi_jp: ", 1, p-j, yi_jp, 1 ); */
	  /*print_matrix( "mu_jp: ", 1, p-j, mu_jp, 1 ); */
	  /*print_matrix( "sigma_jp: ", 1, p-j, sigma_jp, 1 ); */
	  /*Rprintf( "dmvnorm=%f\n", Pijs[i][j][s]); */
	  /*Pij[i][j]+=Pijs[i][j][s];
	  Pjs[j][s]+=Pijs[i][j][s];*/

	  /*free(sigma_jp);
	  free(yi_jp);
	  free(mu_jp);*/
	  free(all_mu);
	  free(mu_p);  
	} 
     }
  }/**/
 /*printing out Pij
 Rprintf("Pij=:");
for (i=0; i<n; i++){
  Rprintf( "i=%d\n", i );
    for (j=0; j<p; j++){
      Rprintf("%2.9f,   ", Pij[i][j]);
    }
   Rprintf( "\n"); 
} */

/*printing out Pjs
Rprintf("Pjs=:");
  for (j=0; j<p; j++){
    Rprintf( "j=%d\n", j );
    for (s=0; s<S; s++){
      Rprintf("%2.9f,   ", Pjs[j][s]);
    }
   Rprintf( "\n"); 
  }*/
 
double *prob_vec, prob_check;
  int *perm, *num_grid;
  prob_vec=(double *) calloc(S, sizeof(double));
  perm=(int*) calloc(S, sizeof(int));
  num_grid=(int*) calloc(1, sizeof(int));
     for (i=0; i<n; i++){
       for (j=0; j<p; j++){	
	 prob_check=0;
	 /*Rprintf("i=%d, j=%d\n", i,j);*/
           for (s=0; s<S; s++){
             if (Pij[i][j]>0) {prob_vec[s]=Pijs[i][j][s]/Pij[i][j];} else {prob_vec[s]=1.0/S;}
             prob_check=prob_check+prob_vec[s];
       }
       /*print_matrix( "prob_vec: ", 1, S, prob_vec, 1 );*/
       if (prob_check<1.0) {for (s=0; s<S; s++) prob_vec[s]=1/S;} 
       
	    int nc = 0;
	    for (s = 0; s < S; s++) if(S * prob_vec[s] > 0.1) nc++;
	    if (nc > 200)
		walker_ProbSampleReplace(S, prob_vec, perm, 1, num_grid);
	    else
		ProbSampleReplace(S, prob_vec, perm, 1, num_grid);      
	    
       eta_mat[i][j]=G[num_grid[0]-1];             
       /* Rprintf("prob_check=%2.9f,  eta_mat[i][j]=%f\n",prob_check, eta_mat[i][j]);*/
      }
   }

   free(prob_vec);
   free(perm);

  prob_vec=(double *) calloc(n, sizeof(double));
  perm=(int*) calloc(n, sizeof(int));
  
/*  for (s=0; s<S; s++){
   for (j=0; j<p; j++){
    for (i=0; i<n; i++){
       if (Pjs[j][s]>0) {prob_vec[i]=Pijs[i][j][s]/Pjs[j][s];} else prob_vec[i]=1.0/n;
   }*/
   /*Rprintf("j=%d,  s=%d",j,s);
       print_matrix( "prob_vec: ", 1, n, prob_vec, 1 );*/   

/*	    int nc = 0;
	    for (i = 0; i < n; i++) if(n * prob_vec[i] > 0.1) nc++;
	    if (nc > 200)
		walker_ProbSampleReplace(n, prob_vec, perm, 1, num_grid);
	    else
		ProbSampleReplace(n, prob_vec, perm, 1, num_grid);      
   

    ystar[S*j+s]=data_mat[num_grid[0]-1][j];
  }
}
   for (j=0; j<p; j++){for (i=0; i<n; i++){ eta[j*n+i]=eta_mat[i][j];}}
*/
  free(prob_vec);
  free(perm);
  free(num_grid);   
  free(eta_i);
}

/* 
 * evaluates and prints the non-zero B-spline basis functions (or their derivatives)
 * at xvals - written for checking  
 */

SEXP splb_check(SEXP Rknots, SEXP Rnknots, SEXP RknotsI, SEXP Rorder, SEXP Reta, SEXP Rm,
                SEXP Rallbeta, SEXP Rn, SEXP Rp, SEXP RS, SEXP RG, SEXP Rall_sigma, SEXP Rdata)
{
    int nknots, order, m, n, i, p, nbeta, S/*, s,s1,j*/ ;
    double *knots, *knotsI, *eta, *all_beta;
    SEXP ans, Reta_mod, Rystar;/**/
    double *ystar, *G, *all_sigma, *data;


    nknots=INTEGER(Rnknots)[0];
    S=INTEGER(RS)[0];

    knots = (double *) calloc(nknots, sizeof(double));
    for(i = 0; i < nknots; i++) knots[i]=REAL(Rknots)[i];
    knotsI = (double *) calloc((nknots-1), sizeof(double));
    for(i = 0; i < (nknots-1); i++) knotsI[i]=REAL(RknotsI)[i];
    G = (double *) calloc(S, sizeof(double));
    for(i = 0; i < S; i++) G[i]=REAL(RG)[i];

    order = INTEGER(Rorder)[0];
    m=INTEGER(Rm)[0];
    n=INTEGER(Rn)[0];
    p=INTEGER(Rp)[0];
    /*Rprintf("order=%d, m=%d, n=%d, p=%d\n", order,m,n,p);*/
    
    double ***Pijs; 
    Pijs=DM3D(n,p,S);
    eta = (double *) calloc(n*p, sizeof(double));
    for(i = 0; i < n*p; i++) eta[i] = REAL(Reta)[i];

    nbeta=p*(p+1)*m/2;
    all_beta=(double *) calloc(nbeta, sizeof(double));
    for(i = 0; i < nbeta; i++) {all_beta[i] = REAL(Rallbeta)[i];}
    all_sigma = (double *) calloc(p, sizeof(double));
    for(i = 0; i < p; i++) all_sigma[i]=REAL(Rall_sigma)[i];
    data = (double *) calloc(n*p, sizeof(double));
    for(i = 0; i < n*p; i++) data[i]=REAL(Rdata)[i];
    ystar = (double *) calloc(S*p, sizeof(double));

    comp_etaY_grid (n, p, S, order, m, nknots, knots, knotsI,  G, all_beta, all_sigma, data, eta, ystar, Pijs);/**/

  /*printing out Pij
  Rprintf("Pijs=:");
  for(s=0; s<S; s++){
    Rprintf("s=%d\n", s);
    for (i=0; i<n; i++){
  for (j=0; j<p; j++){
	  Rprintf("Pijs[%d][%d][%d]=%2.9f,   ", i,j,s, Pijs[i][j][s]);
	}
      Rprintf( "\n"); 
  } } */


    PROTECT(ans = allocVector(VECSXP, 2));
    PROTECT(Reta_mod = allocMatrix(REALSXP, n, p));
    PROTECT(Rystar = allocMatrix(REALSXP, S, p));
    for(i = 0; i < n*p; i++) {REAL(Reta_mod)[i] = eta[i];}
    for(i = 0; i < S*p; i++) {REAL(Rystar)[i] = ystar[i];}
    SET_VECTOR_ELT(ans, 0, Reta_mod);
    SET_VECTOR_ELT(ans, 1, Rystar);
    UNPROTECT(3);

    
    
   /*PROTECT(ans = allocVector(VECSXP, S));
    for (s1=0; s1<S; s1++){
      PROTECT(Rystar = allocMatrix(REALSXP, n, p));
      for(j = 0; j < p; j++) { for (i=0; i<n; i++) {
        REAL(Rystar)[j*n+i]=Pijs[i][j][s1];
        SET_VECTOR_ELT(ans, s1, Rystar);}
      }
      UNPROTECT(1);
    }
    UNPROTECT(1); */


    free(eta);
    free(all_beta);
    free(all_sigma);
    free(G);
    free(data);
    free(knotsI);
    free(knots);
    free(ystar);
    Dfree3d(Pijs);
    return(ans);
}

/*
 * 
 */
void mu_B_fct (int n, int p, int order, int m, int nknots, double *knots, double *knotsI,  
                double *all_beta, double *eta, double *Ball_mat, double *mu)
{

  int i, j,k, supi,l;
  double *B, *BI, *Ball, *all_mu, *mu_p, *etai; 
  int *offsets, *offsetsI, *derivs;

for (i=0; i<n; i++){
  Ball=(double*)calloc(m*p, sizeof(double));
  all_mu=(double *) calloc(p*(p+1)/2, sizeof(double));
  B=(double*)calloc(order*p, sizeof(double));
  BI=(double*)calloc(order*p, sizeof(double));
  mu_p=(double *) calloc(p, sizeof(double));
  etai=(double *) calloc(p, sizeof(double));
  derivs=(int*)calloc(p, sizeof(int));    
  offsets=(int*)calloc(p, sizeof(int));
  offsetsI=(int*)calloc(p, sizeof(int));
  
  for (k=0; k<p; k++) etai[k]=eta[k*n+i];
  spline_basis(knotsI, (nknots-1), order, etai, p, derivs, p, offsetsI, BI);
  spline_basis(knots, nknots, order, etai, p, derivs, p, offsets, B);
  for (l=0; l<order; l++){
      Ball[offsetsI[0]+l]=BI[l];
    }  
  for (j=1; j<p; j++){
    if (offsets[j]!=0) {
      for (l=0; l<order; l++){Ball[m*j+offsets[j]+l-1]=B[j*order+l];}
    } else {for (l=1; l<order; l++){ Ball[m*j+offsets[j]+l-1]=B[j*order+l];}}
  } 
  for (j=0;j<p; j++){
    for (k=0; k<=j; k++){  
      for (supi=k*m; supi<(k+1)*m; supi++) all_mu[j*(j+1)/2+k]+=Ball[supi] * all_beta[j*(j+1)*m/2+supi];
    }
  }
  for (j=0; j<p; j++){for (k=0;k<j+1; k++){ mu_p[j]=mu_p[j]+all_mu[j*(j+1)/2+k]; }}
  for (k=0; k<p; k++) mu[k*n+i]=mu_p[k];
  for (k=0; k<m*p; k++) Ball_mat[k*n+i]=Ball[k];
  free(Ball);
  free(all_mu);
  free(mu_p);
  free(B);
  free(BI);
  free(etai);
  free(offsets);
  free(offsetsI);
  free(derivs);
}
}


/* this function checks whether we compute correctly parameters for betaj posterior distribution*/

SEXP betaj_post_check(SEXP Rknots, SEXP Rnknots, SEXP RknotsI, SEXP Rorder, SEXP Reta, SEXP Rm,
                SEXP Rallbeta, SEXP Rn, SEXP Rp, SEXP RS, SEXP RG, SEXP Rall_sigma, SEXP Rdata, SEXP Rtau)
{

    int nknots, order, m, n, i, j, p, nbeta, S, mj1/*, s,s1,j*/ ;
    double *knots, *knotsI, *eta, *all_beta;
    SEXP ans, Rcov, Rhalfcov, Rmean;/**/
    double *all_sigma, *data, *Ball_mat, *Bj_mat, *tauj;
    double sigj,*mu, *meanj, *covj, *half_covj, *yj; 


    nknots=INTEGER(Rnknots)[0];
    S=INTEGER(RS)[0];

    knots = (double *) calloc(nknots, sizeof(double));
    for(i = 0; i < nknots; i++) knots[i]=REAL(Rknots)[i];
    knotsI = (double *) calloc((nknots-1), sizeof(double));
    for(i = 0; i < (nknots-1); i++) knotsI[i]=REAL(RknotsI)[i];

    order = INTEGER(Rorder)[0];
    m=INTEGER(Rm)[0];
    n=INTEGER(Rn)[0];
    p=INTEGER(Rp)[0];
    eta = (double *) calloc(n*p, sizeof(double));
    for(i = 0; i < n*p; i++) eta[i] = REAL(Reta)[i];
    nbeta=p*(p+1)*m/2;
    all_beta=(double *) calloc(nbeta, sizeof(double));
    for(i = 0; i < nbeta; i++) {all_beta[i] = REAL(Rallbeta)[i];}
    all_sigma = (double *) calloc(p, sizeof(double));
    for(i = 0; i < p; i++) all_sigma[i]=REAL(Rall_sigma)[i];
    
    data = (double *) calloc(n*p, sizeof(double));
    for(i = 0; i < n*p; i++) data[i]=REAL(Rdata)[i];
    
    Ball_mat=(double*)calloc(m*p*n, sizeof(double));
    mu=(double*)calloc(p*n, sizeof(double));
        
    mu_B_fct(n, p, order,m, nknots,knots, knotsI,  all_beta, eta, Ball_mat, mu);
    /*print_matrix("Ball_mat", n,m*p, Ball_mat, m*p);*/

    PROTECT(ans = allocVector(VECSXP, 9));
    
    for (j=0; j<p; j++) {
      mj1=m*(j+1);
      Bj_mat=(double*)calloc(mj1*n, sizeof(double));
      for (i=0; i<mj1*n; i++) Bj_mat[i]=Ball_mat[i];
      
      tauj=(double*)calloc(j, sizeof(double));
      for (i=0; i<(j+1); i++) tauj[i]=REAL(Rtau)[i];
      /*Rprintf("j+1=%d", j+1);
      print_matrix("tauj", (j+1), 1,tauj, 1);*/
      
      yj=(double*)calloc(n, sizeof(double));
      for (i=0; i<n; i++) yj[i]=data[n*j+i];
      /*print_matrix("yj", n, 1,yj, 1);*/
      
      meanj=(double*)calloc(mj1, sizeof(double));
      covj=(double*)calloc(mj1*mj1, sizeof(double));
      half_covj=(double*)calloc(mj1*mj1, sizeof(double));
      
      sigj=all_sigma[j];
      post_betaj_mean_cov(Bj_mat, tauj, yj,  sigj,  j,  m,  n, meanj, covj, half_covj);
      
      PROTECT(Rcov = allocMatrix(REALSXP, mj1, mj1));
      PROTECT(Rhalfcov = allocMatrix(REALSXP, mj1, mj1));
      PROTECT(Rmean = allocVector(REALSXP, mj1));
      for(i = 0; i < mj1*mj1; i++) {REAL(Rcov)[i] = covj[i]; REAL(Rhalfcov)[i] = half_covj[i];}
      for(i = 0; i < mj1; i++) {REAL(Rmean)[i] = meanj[i];}
      SET_VECTOR_ELT(ans, j*3, Rcov);
      SET_VECTOR_ELT(ans, j*3+1, Rhalfcov);
      SET_VECTOR_ELT(ans, j*3+2, Rmean);      
      free(Bj_mat);
      free(tauj);
      free(yj);
      free(meanj);
      free(covj);
      free(half_covj);
      UNPROTECT(3);
    }
UNPROTECT(1);   

    
    free(eta);
    free(all_beta);
    free(all_sigma);
    free(data);
    free(knotsI);
    free(knots);
    free(Ball_mat);
    free(mu);
    return(ans);
}

/*******************************************************************************************************/

SEXP gibbs_BALVM_naive(SEXP Rtot_draws, SEXP Rknots, SEXP RNknots, SEXP RknotsI, SEXP Rorder,
          SEXP Rn, SEXP Rm, SEXP Rp, SEXP RGrid_length, SEXP Rgrid, SEXP Rdata, SEXP Rab_sigma, SEXP Rab_tau, SEXP Reta_start, SEXP Rallbeta_start)
{
 int i, k, j, s, m, p, n, order, grid_length;
 int pn, mj1, t, tot_draws, Nknots;
 double *data, *ystar, *knots, *knotsI;
 double *eta, *etastar, *all_beta, *betaj;
 double *mu, *Ball_mat;
 double a_sigma, b_sigma, a_tau, b_tau;
 double sig_j, tau_kj;
 double *sigma_vec, *tauj, *yj;
 double *Bj_mat, *muj, *covj, *half_covj, *betakj, *mu_betaj;
 double *grid;
 SEXP ans;
 
 Nknots=INTEGER(RNknots)[0];
 m=INTEGER(Rm)[0];
 p=INTEGER(Rp)[0];
 order=INTEGER(Rorder)[0];
 n=INTEGER(Rn)[0];
 pn=p*n;
 grid_length=INTEGER(RGrid_length)[0];
 tot_draws=INTEGER(Rtot_draws)[0];
 a_sigma=REAL(Rab_sigma)[0];
 b_sigma=REAL(Rab_sigma)[1];
 a_tau=REAL(Rab_tau)[0];
 b_tau=REAL(Rab_tau)[1];
 
 data=(double *) calloc(pn, sizeof(double));
 eta=(double *) calloc(pn, sizeof(double));
 for (i=0; i<pn; i++) {data[i]=REAL(Rdata)[i]; eta[i]=REAL(Reta_start)[i];}
 
 grid=(double *) calloc(grid_length, sizeof(double));
 for (i=0; i<grid_length; i++) grid[i]=REAL(Rgrid)[i];
 
 all_beta=(double *) calloc(p*(p+1)*m/2, sizeof(double));
 for (i=0; i<p*(p+1)*m/2; i++) all_beta[i]=REAL(Rallbeta_start)[i];
 
 knots=(double *) calloc(Nknots, sizeof(double));
 for (i=0; i<Nknots; i++) knots[i]=REAL(Rknots)[i];
 knotsI=(double *) calloc(Nknots-1, sizeof(double));
 for (i=0; i<(Nknots-1); i++) knotsI[i]=REAL(RknotsI)[i];
 
 mu=(double *) calloc(pn, sizeof(double));
 Ball_mat=(double *) calloc(n*m*p, sizeof(double));
 mu_B_fct (n, p, order, m, Nknots, knots, knotsI,  all_beta, eta, Ball_mat, mu);
 /*print_matrix( "mu matrix should be n x p", n, p, mu, n);*/
 sigma_vec=(double *) calloc(p, sizeof(double));
 
 double ***Pijs; 
 Pijs=DM3D(n,p,grid_length);


GetRNGstate();
 /*starting draws*/
for (t=0; t<tot_draws; t++)
 {
  for (j=0; j<p; j++)
  {
    mj1=m*(j+1);
    Bj_mat=(double*)calloc(mj1*n, sizeof(double));
    for (i=0; i<mj1*n; i++) Bj_mat[i]=Ball_mat[i];

    tauj=(double*)calloc(j, sizeof(double));
    yj=(double*)calloc(n, sizeof(double));
    for (i=0; i<n; i++) yj[i]=data[n*j+i];

    muj=(double*)calloc(n, sizeof(double));
    for (i=0; i<n; i++) muj[i]=mu[j*n+i];
    /*printf("j=%d  j*n=%d  \n", j, j*n);
    print_matrix( "muj", 1, n, muj, 1);*/

    betaj=(double*)calloc(mj1, sizeof(double));
    for (s=0; s<mj1; s++) betaj[s]=all_beta[j*(j+1)*m/2+s];

    sig_j=sim_sigma_j(yj, muj, n, a_sigma, b_sigma);
    sigma_vec[j]=sig_j;
    for (k=0; k<j; k++)
    {
      betakj=(double*)calloc(m, sizeof(double));
      for (s=0; s<m; s++) betakj[s]=betaj[k*m+s];
      tau_kj=sim_taukj_naive(a_tau, b_tau, betakj, m);
      tauj[k]=tau_kj;
      free(betakj);
    }
    
    mu_betaj=(double*)calloc(mj1, sizeof(double));
    covj=(double*)calloc(mj1*mj1, sizeof(double));
    half_covj=(double*)calloc(mj1*mj1, sizeof(double));   
    post_betaj_mean_cov(Bj_mat, tauj, yj, sig_j, j, m, n, muj, covj, half_covj);
    sim_betaj(j, m, mu_betaj, half_covj, betaj);
    for (s=0; s<mj1; s++) all_beta[j*(j+1)*m/2+s]=betaj[s];/**/
    
    free(Bj_mat);
    free(tauj);
    free(yj);
    free(muj);
    free(covj);
    free(half_covj);
    free(betaj);
    free(mu_betaj); /**/
  }  /* end j loop */

 comp_etaY_grid (n, p, grid_length, order, m, Nknots, knots, knotsI, 
                grid, all_beta, sigma_vec, data, eta, ystar, Pijs);

} 
/* end draws*/


PutRNGstate(); 
 free(data);
 free(eta);
 free(grid);
 free(all_beta);
 free(knots);
 free(knotsI);
 free(mu);
 free(Ball_mat);
 free(sigma_vec);
 Dfree3d(Pijs);
 PROTECT(ans = allocVector(REALSXP,1));
 UNPROTECT(1);
 return ans;
}

/******************************************************************************************************/

/* evaluate the non-zero B-spline basis functions (or their derivatives)
 * at xvals.  */
 
SEXP B_mu_check(SEXP Rknots, SEXP Rnknots, SEXP RknotsI, SEXP Rorder, SEXP Reta, SEXP Rm,
                SEXP Rallbeta, SEXP Rn, SEXP Rp, SEXP RS, SEXP RG, SEXP Rall_sigma, SEXP Rdata)
{

    int nknots, order, m, n, i, p, nbeta, S/*, s,s1,j*/ ;
    double *knots, *knotsI, *eta, *all_beta;
    SEXP ans, RBall_mat, Rmu;
    double *all_sigma, *data, *Ball_mat, *mu;


    nknots=INTEGER(Rnknots)[0];
    S=INTEGER(RS)[0];

    knots = (double *) calloc(nknots, sizeof(double));
    for(i = 0; i < nknots; i++) knots[i]=REAL(Rknots)[i];
    knotsI = (double *) calloc((nknots-1), sizeof(double));
    for(i = 0; i < (nknots-1); i++) knotsI[i]=REAL(RknotsI)[i];

    order = INTEGER(Rorder)[0];
    m=INTEGER(Rm)[0];
    n=INTEGER(Rn)[0];
    p=INTEGER(Rp)[0];
    eta = (double *) calloc(n*p, sizeof(double));
    for(i = 0; i < n*p; i++) eta[i] = REAL(Reta)[i];
    nbeta=p*(p+1)*m/2;
    all_beta=(double *) calloc(nbeta, sizeof(double));
    for(i = 0; i < nbeta; i++) {all_beta[i] = REAL(Rallbeta)[i];}
    all_sigma = (double *) calloc(p, sizeof(double));
    for(i = 0; i < p; i++) all_sigma[i]=REAL(Rall_sigma)[i];
    data = (double *) calloc(n*p, sizeof(double));
    for(i = 0; i < n*p; i++) data[i]=REAL(Rdata)[i];
    Ball_mat=(double*)calloc(m*p*n, sizeof(double));
    mu=(double*)calloc(p*n, sizeof(double));
    
    mu_B_fct(n, p, order,m, nknots,knots, knotsI,  all_beta, eta, Ball_mat, mu);

    PROTECT(ans = allocVector(VECSXP, 2));
    PROTECT(RBall_mat = allocMatrix(REALSXP, n, m*p));
    PROTECT(Rmu = allocMatrix(REALSXP, n, p));
    for(i = 0; i < n*m*p; i++) {REAL(RBall_mat)[i] = Ball_mat[i];}
    for(i = 0; i < n*p; i++) {REAL(Rmu)[i] = mu[i];}
    SET_VECTOR_ELT(ans, 0, RBall_mat);
    SET_VECTOR_ELT(ans, 1, Rmu);
    UNPROTECT(3);

    
    free(eta);
    free(all_beta);
    free(all_sigma);
    free(data);
    free(knotsI);
    free(knots);
    free(Ball_mat);
    free(mu);
    return(ans);
}

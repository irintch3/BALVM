#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#define LOW -1.0e15

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("splines", String)
#else
#define _(String) (String)
#endif

#define HIGH 1.0e15

/********************************/

long fsafewrite(ptr,size,n,stream)
     double *ptr;size_t size;long n;FILE *stream;
 /*
 * safely write into a file
 */
{ long i,j,k=0L;
  for (i=0;i<(n/32000L);i++)
  k+=fwrite(ptr+i*32000L,size,(size_t)32000L,stream);
  j=n%32000L;
  k+=fwrite(ptr+i*32000L,size,(size_t)j,stream);
  return(k);
}

long fsaferead(ptr,size,n,stream)
     double *ptr;size_t size;long n;FILE *stream;
/*
 * safely read from a file
 */
{ long i,j=32000L,k=0L;
  for (i=0;i<(n/32000L);i++)
  k+=fread(ptr+i*32000L,size,(size_t)j,stream);
  j=n%32000L;
  k+=fread(ptr+i*32000L,size,(size_t)j,stream);
  return(k);
}
/******************************************************************************/

void computeWMat(double *eta, int p, int order, int m, int nknots, double *knots, double *knotsI, double *all_beta, double *all_mu)
/*all_mu should be an array of zeros*/
{
/*
m is the degree of freedom of splines and number of components in beta_{kj}
knots and knotsI is for splines without intercept (but the first coordinates is the intercept function and should be excluded)
number of knots in knotsI is knots-1
order is the order of the spline which is the degree+1
n is the number od observations
*/

/*******************************************************************/

void all_mu_comp(double *eta, double *G, double *all_beta, double *all_mu)
/*all_mu should be an array of zeros*/
{
/*
m is the degree of freedom of splines and number of components in beta_{kj}
knots and knotsI is for splines without intercept (but the first coordinates is the intercept function and should be excluded)
number of knots in knotsI is knots-1
order is the order of the spline which is the degree+1
n is the number od observations
*/

  int i, j,k, supi;
  double *B, *BI, *Ball;
  int *offsets, *offsetsI, *derivs;

  B=(double*)calloc(order*p, sizeof(double));
  BI=(double*)calloc(order*p, sizeof(double));
  Ball=(double*)calloc(m*p, sizeof(double));

  derivs=(int*)calloc(p, sizeof(int));
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

/************************************************************************/
/*
 *  Unequal Probability Sampling.
 *
 *  Modelled after Fortran code provided by:
 *    E. S. Venkatraman <venkat@biosta.mskcc.org>
 *  but with significant modifications in the
 *  "with replacement" case.
 */

/* Unequal probability sampling; with-replacement case */

void ProbSampleReplace(int n, double *p, int *perm, int nans, int *ans)
{
    double rU;
    int i, j;
    int nm1 = n - 1;

    /* record element identities */
    for (i = 0; i < n; i++)
  	perm[i] = i + 1;

    /* sort the probabilities into descending order */
    revsort(p, perm, n);

    /* compute cumulative probabilities */
    for (i = 1 ; i < n; i++)
  	p[i] += p[i - 1];

    /* compute the sample */
    for (i = 0; i < nans; i++) {
  	rU = unif_rand();
	   for (j = 0; j < nm1; j++) {
  	    if (rU <= p[j])
 		break;
	}
	ans[i] = perm[j];
    }
}


/* use example to run in R under linux

    dyn.load("fcts_Eta_B_win.so")
    Rtot_draws <- as.integer(18)#0000
    Rdens=0
    Rdrws.step <- 1

    mu.all.my <- .Call("gibbs_BALVM_naive", Rtot_draws, Aknots, as.integer(nknots), AknotsI, as.integer(ord),
                       as.integer(n), as.integer(m), as.integer(p), as.integer(S), grid.eta, data, c(1,0.05),
                       c(1,0.05), eta, Rallbeta, as.integer(Rdens),as.integer(Rdrws.step) )

    dyn.unload("fcts_Eta_B_win.so")
*/

double dmvnorm(int dim, double *all_sigma, double *y, double *mu)
{
  /*
  compute multivarite normal density
  all_sigma contains the diag vector of the covariance matrix
  */
  int i;
  double val;
  double sum_sq, lndet;
  sum_sq=0;
  lndet=0;

  for (i=0; i<dim; i++){sum_sq=sum_sq+pow(y[i]-mu[i], 2.0)/(2*all_sigma[i]); lndet=lndet+0.5*log(all_sigma[i]*2*PI);}
  val=exp(-sum_sq-lndet);
  return val;
}



/******************************************************************************/
double ***DM3D(int n,int row,int col) {
/*
* this function allocates memory to 3-dimensional array
*/
  int i,j;
  /*int ibase=0;*/
  double ***mem;
  double **prow;
  double *pdata;

  pdata= (double *) calloc(n*row*col, sizeof(double));
  if ( pdata==(double *) NULL ) {
     fprintf(stderr, "No heap space for 3D data\n");
     exit(1);
  }

  prow=(double **) calloc(n*row, sizeof(double *));
  if ( prow==(double **) NULL ) {
     fprintf(stderr, "No heap space for 3D data\n");
     exit(1);
  }
  mem=(double ***) calloc(n,sizeof(double **));
  if ( mem==(double ***) NULL ) {
     fprintf(stderr, "No heap space for 3D data\n");
     exit(1);
  }

  for (i=0;i<n;i++) {
    for ( j=0;j<row; j++) {
      prow[j]=pdata;
      pdata+=col;
    }
    mem[i]=prow;
    prow+=row;
  }
  return mem;
}
/******************************************************************************/

void Dfree3d(double ***pa) {
/*
* this function frees memory taken by 3-dimensional array
*/

/* free the data */
free(**pa);
/* free the row pointers */
free(*pa);
/* free the 3D pointers */
free(pa);
}

/******************************************************************************/
double **DM2D(int n,int row) {
/*
* this function allocates memory to a 2-dimensional array
*/
  int i;
  /*int ibase=0;*/
  double **mem;
  double *pdata;

  pdata= (double *) calloc(n*row, sizeof(double));
  if ( pdata==(double *) NULL ) {
     fprintf(stderr, "No heap space for 2D data\n");
     exit(1);
  }

  mem=(double **) calloc(n,sizeof(double *));
  if ( mem==(double **) NULL ) {
     fprintf(stderr, "No heap space for 2D data\n");
     exit(1);
  }

  for (i=0;i<n;i++) {
    mem[i]=pdata;
    pdata+=row;
  }
  return mem;
}
/******************************************************************************/

void Dfree2d(double **pa) {
/*
* this function frees memory taken by 2-dimensional array
*/
free(*pa);
free(pa);
}

/******************************************************************************/

void comp_etaY_grid (int n, int p, int S, int order, int m, int Nknots,
		     double *knots, double *knotsI,  double *G, double *all_beta, double *all_sigma,
		     double *data, double *eta, double *ystar, double ***Pijs, int dens_comp)
{
/*
* This procedure computes and simulates the new eta and ystar
* in the paper eta is denoted by u
*/
  double *all_mu, *yi_jp,  *mu_p, *mu_jp, *eta_i, *sigma_jp;
  double **Pij, **Pjs;
  int s, i,j,k, j1;
  double *prob_vec, prob_check;
  int *perm, *num_grid;

  Pij=DM2D(n,p);
  Pjs=DM2D(p,S);
  for (j=0; j<p;j++){
    for (s=0; s<S; s++) Pjs[j][s]=0;
    for (i=0; i<n; i++) Pij[i][j]=0;
  }

  eta_i=(double *) malloc(p*sizeof(double));
  all_mu=(double *) calloc(p*(p+1), sizeof(double)/2);
  mu_p=(double *) calloc(p, sizeof(double));


for (i=0; i<n; i++){
  for (j=0; j<p; j++){
    for (s=0; s<S; s++){
  	  for (k=0; k<p; k++) eta_i[k]=eta[k*n+i];
  	  eta_i[j]=G[s];
  	  all_mu_comp(eta_i, p, order, m, Nknots, knots, knotsI, all_beta, all_mu);
  	  for (j1=0; j1<p; j1++){for (k=0;k<j1+1; k++){ mu_p[j1]=mu_p[j1]+all_mu[j1*(j1+1)/2+k]; all_mu[j1*(j1+1)/2+k]=0;}}
  	  sigma_jp=(double *) calloc((p-j),sizeof(double));
  	  yi_jp=(double *) calloc((p-j),sizeof(double));
  	  mu_jp=(double *) calloc((p-j),sizeof(double));
  	  for (k=0; k<p-j; k++) { yi_jp[k]=data[(j+k)*n+i]; sigma_jp[k]=all_sigma[j+k]; mu_jp[k]=mu_p[j+k];}
  	  Pijs[i][j][s]=dmvnorm(p-j, sigma_jp, yi_jp, mu_jp);
      memset(mu_p, 0, p*sizeof(double));
  	  Pij[i][j]+=Pijs[i][j][s];
  	  if (dens_comp==1) Pjs[j][s]+=Pijs[i][j][s];
  	  free(sigma_jp);
  	  free(yi_jp);
  	  free(mu_jp);
	   }
    }
  }
free(all_mu);
free(mu_p);
free(eta_i);

/*printing out Pij - used in debugging
Rprintf("Pij=:");
for (i=0; i<n; i++){
  Rprintf( "i=%d\n", i );
    for (j=0; j<p; j++){
      Rprintf("%2.9f,   ", Pij[i][j]);
    }
   Rprintf( "\n");
}  */

/* printing out Pjs  - used in debugging
Rprintf("Pjs=:");
  for (j=0; j<p; j++){
    Rprintf( "j=%d\n", j );
    for (s=0; s<S; s++){
      Rprintf("%2.9f,   ", Pjs[j][s]);
    }
   Rprintf( "\n");
  }*/


  perm=(int*) calloc(S, sizeof(int));
  num_grid=(int*) calloc(1, sizeof(int));
  prob_vec=(double *) malloc(S*sizeof(double));

     for (i=0; i<n; i++){
       for (j=0; j<p; j++){
           prob_check=0;
           if (Pij[i][j]!=0)
           {
            for (s=0; s<S; s++)
              {
                prob_vec[s]=Pijs[i][j][s]/Pij[i][j];
                /*Rprintf("Pijs[%d][%d][%d]=%f, Pij[%d][%d]=%f, Pijs[%d][%d][%d]/Pij[%d][%d]=%f,  prob_vec[%d]=%f\n",i,j,s,Pijs[i][j][s], i,j, Pij[i][j],  i,j,s,i,j,Pijs[i][j][s]/Pij[i][j], s, prob_vec[s]);*/
                prob_check=prob_check+prob_vec[s];
              }
           } else {for (s=0; s<S; s++) prob_vec[s]=1.0/S;}
          /*print_matrix(" comp_etaY_grid prob_vec=", S,1,prob_vec,S );*/
	  ProbSampleReplace(S, prob_vec, perm, 1, num_grid);
	  /*Rprintf("num_grid[0]-1=%d\n", num_grid[0]-1);*/
          eta[i+j*n]=G[num_grid[0]-1];
          }
   }
  free(prob_vec);
  free(perm);

if (dens_comp==1)
{
  prob_vec=(double *) malloc(n*sizeof(double));
  perm=(int*) calloc(n, sizeof(int));

  for (s=0; s<S; s++){
   for (j=0; j<p; j++){
    if (Pjs[j][s]!=0) { for (i=0; i<n; i++) prob_vec[i]=Pijs[i][j][s]/Pjs[j][s];} else prob_vec[i]=1.0/n;
    /*Rprintf("num_grid[0]-1=%d\n", num_grid[0]-1);*/
    ProbSampleReplace(n, prob_vec, perm, 1, num_grid);
    ystar[S*j+s]=data[num_grid[0]-1+j*n];
  }
 }
  free(prob_vec);
  free(perm);
}

  free(num_grid);
  free(*Pij);
  free(Pij);
  free(*Pjs);
  free(Pjs);

}


 /******************************************************************************/
/*
* this function is taken from Simon WOod package mgcv
*/

void getXtX(double *X,int *r,int *c, double *XtX)
/* form X'X (nearly) as efficiently as possible
*   r is number of rows,
*   c is number of columns
*/
{ double *p0,*p1,*p2,*p3,*p4,x;
  int i,j;
  for (p0=X,i=0;i<*c;i++,p0 += *r)
  for (p1=X,j=0;j<=i;j++,p1 += *r) {
    for (x=0.0,p2=p0,p3=p1,p4=p0 + *r;p2<p4;p2++,p3++) x += *p2 * *p3;
    XtX[i + j * *c] = XtX[j + i * *c] = x;
  }
}

/***************************************************************************/

void transpose_mat(double *mat, int ncol, int nrow, double *resmat)
{
 int j=0, i=0;
       for(i = 0; i < nrow; i++)    /* column number for the element of resmat */
         for (j = 0; j < ncol; j++)     /* line number for the element of resmat */
           {
              resmat[i* ncol+j]=mat[j* nrow+i];
           }
}

/***************************************************************************/

void matprod(double *mat1, int ncol1, int nrow1, double *mat2, int ncol2, int nrow2, double *resmat)
/*
computes product of two matrices ncol1 should be equal nrow2
otherwise nothing happens
*/
 {
int i=0, j=0, k=0;
double sum;
 if (ncol1==nrow2) {
     for (j = 0; j < ncol2; j++)     /* column number for the element of resmat */
       for(i = 0; i < nrow1; i++)    /* line number for the element of resmat */
          {
             sum=0;
             for (k = 0; k < ncol1; k++)  sum=sum+ mat1[k* nrow1+i] * mat2[j* nrow2+k];
              resmat[j* nrow1+i]=sum;
          }
    }
}

/**************************************************************************************************************************************************************/

/*********************************************************************************/


void mu_B_fct (int n, int p, int order, int m, int nknots, double *knots,
	       double *knotsI,  double *all_beta, double *eta, double *Ball_mat,
	       double *mu, FILE *out_allmu, int save)
{
/*
* this procedure computes the overall matrix B (Ball_mat) and the total values of mu_1, ..., mu_p (mu)
*/

  int i, j,k, supi,l;
  double *B, *BI, *Ball, *all_mu, *mu_p, *etai;
  int *offsets, *offsetsI, *derivs;
  Ball=(double*)calloc(m*p, sizeof(double));
  all_mu=(double *) calloc(p*(p+1)/2, sizeof(double));
  B=(double*)calloc(order*p, sizeof(double));
  BI=(double*)calloc(order*p, sizeof(double));
  mu_p=(double *) calloc(p, sizeof(double));
  etai=(double *) calloc(p, sizeof(double));
  derivs=(int*)calloc(p, sizeof(int));
  offsets=(int*)calloc(p, sizeof(int));
  offsetsI=(int*)calloc(p, sizeof(int));

  for (i=0; i<n; i++)
  {
    memset( Ball, 0, m*p*sizeof(double));
    memset( all_mu,0, p*(p+1)*sizeof(double)/2);
    memset( B,0, order*p*sizeof(double));
    memset( BI,0, order*p*sizeof(double));
    memset( mu_p,0, p*sizeof(double));
    memset( derivs,0, p*sizeof(int));
    memset( offsets,0, p*sizeof(int));
    memset( offsetsI,0, p*sizeof(int));
    for (k=0; k<p; k++) etai[k]=eta[k*n+i];

    spline_basis(knotsI, (nknots-1), order, etai, p, derivs, p, offsetsI, BI);
    spline_basis(knots, nknots, order, etai, p, derivs, p, offsets, B);

    for (l=0; l<order; l++) Ball[offsetsI[0]+l]=BI[l];

    for (j=1; j<p; j++)
    {
      if (offsets[j]!=0)
      {
	     for (l=0; l<order; l++)Ball[m*j+offsets[j]+l-1]=B[j*order+l];
      }
      else {for (l=1; l<order; l++) Ball[m*j+offsets[j]+l-1]=B[j*order+l];}
    }

    for (j=0;j<p; j++)
    {
      for (k=0; k<=j; k++)
      {
      	for (supi=k*m; supi<(k+1)*m; supi++) all_mu[j*(j+1)/2+k]+=Ball[supi] * all_beta[j*(j+1)*m/2+supi];
      }
    }

    for (j=0; j<p; j++){for (k=0;k<j+1; k++){ mu_p[j]=mu_p[j]+all_mu[j*(j+1)/2+k]; }}
    for (k=0; k<p; k++) mu[k*n+i]=mu_p[k];
    for (k=0; k<m*p; k++) Ball_mat[k*n+i]=Ball[k];
    if (save==1) fsafewrite(all_mu,sizeof(double),(p+1)*p/2,  out_allmu); /**/
}
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


/*********************************************************************************/

/*******************************************************************************************************/


void one_step_gibbs_naive(int m, int n, int p, double *Ball_mat, double *data, double *mu, double a_sigma,
			  double b_sigma, double a_tau, double b_tau,double *all_beta, double *sigma_vec )
{
/*
* one iteration of Gibbs sampler for BALVM with naive prior
*/
 int j,i,s, mj1, k;
 double sig_j, tau_kj;
 double *Bj_mat, *betaj,*betakj, *tauj, *mu_betaj, *half_covj, *yj, *muj;
 double scl, shp;

 yj=(double*)malloc(n*sizeof(double));
 muj=(double*)malloc(n*sizeof(double));
 betakj=(double*)malloc(m*sizeof(double));

 for (j=0; j<p; j++)
  {
    mj1=m*(j+1);
    
    Bj_mat=(double*)calloc(mj1*n, sizeof(double));
    betaj=(double*)calloc(mj1, sizeof(double));
    tauj=(double*)calloc((j+1), sizeof(double));
    mu_betaj=(double*)calloc(mj1, sizeof(double));
    half_covj=(double*)calloc(mj1*mj1, sizeof(double));
    
    for (i=0; i<mj1*n; i++) Bj_mat[i]=Ball_mat[i];
    for (i=0; i<n; i++) yj[i]=data[n*j+i];
    for (i=0; i<n; i++) muj[i]=mu[j*n+i];
    for (s=0; s<mj1; s++) betaj[s]=all_beta[j*(j+1)*m/2+s];

    /*sig_j simulation*/
    shp=a_sigma+0.5*n;
    scl=b_sigma;  for (i=0; i<n; i++) scl=scl+0.5*pow(yj[i]-muj[i],2);
    sig_j=0.1;/*1/rgamma(shp, 1/scl);*/
    sigma_vec[j]=sig_j;
    /*Rprintf("j=%d, sig_j=%f\n", j, sig_j);
    end sig_j simulation*/

    /*tau_kj simulation*/
    shp=a_tau+0.5*m;
    for (k=0; k<j+1; k++)
    {
      for (s=0; s<m; s++) betakj[s]=betaj[k*m+s];
      scl=b_tau; for (i=0; i<m; i++) scl=scl+0.5*pow(betakj[i],2);
      tau_kj=1/rgamma(shp, 1/scl);
      tauj[k]=tau_kj;
    }
    /*end tau_kj simulation*/
    /*posterior computation for betaj*/
      double *BtB, *Bt, *sigtauj1, *diagcovj, *covj;
      int get_inv;

      BtB=(double *) calloc(mj1*mj1, sizeof(double));
      Bt=(double *) calloc(mj1*n, sizeof(double));
      sigtauj1=(double *) calloc((j+1), sizeof(double));
      diagcovj=(double *) calloc(mj1,sizeof(double));
      covj=(double*)calloc(mj1*mj1, sizeof(double));

      getXtX(Bj_mat,&n,&mj1, BtB);
      for (i=0; i<(j+1); i++) sigtauj1[i]=sig_j/tauj[i];
      for (i=0; i<(j+1); i++) {for (k=0; k<m;k++) BtB[mj1*m*i+mj1*k+m*i+k]=BtB[mj1*m*i+mj1*k+m*i+k]+sigtauj1[i];}
      get_inv=1;
      qr_ldet_inv(BtB,&mj1,covj,&get_inv, diagcovj);
      for (i=0; i<mj1*mj1; i++) half_covj[i]=covj[i];  /*half_covj=covj*/
      free(sigtauj1);
      free(diagcovj);
      free(BtB);

      BtB=(double *) calloc(mj1*n, sizeof(double));
      transpose_mat(Bj_mat, mj1, n, Bt);
      matprod(covj, mj1, mj1, Bt,  n, mj1, BtB); /*BtB=(BtB)^{-1}%*%t(Bj_mat)*/
      matprod(BtB,  n, mj1, yj, 1, n, mu_betaj); /*mu_betaj=(BtB)^{-1}%*%t(Bj_mat)%*%yj */
      free(covj);
      free(BtB);
      free(Bt);

      mroot(half_covj,&mj1,&mj1);/*half_covj=t((covj)^0.5)*/
      BtB=(double *) calloc(mj1*mj1, sizeof(double));
      transpose_mat(half_covj, mj1, mj1, BtB); /*half_covj=(covj)^0.5*/
      for (i=0; i<mj1*mj1; i++) half_covj[i]=BtB[i];
      free(BtB);
      /*end posterior computation for betaj*/
      /*simulation of betaj*/
      double *hlp;
      hlp=(double *) calloc(mj1, sizeof(double));
      for (k=0; k<mj1; k++) betaj[k]=rnorm(0, 1);
      matprod(half_covj, mj1, mj1, betaj, 1, mj1, hlp);
      for (k=0; k<mj1; k++) betaj[k]=hlp[k]+mu_betaj[k];
      free(hlp);
      free(mu_betaj);
      free(half_covj);
      /*end simulation of betaj*/

      for (s=0; s<mj1; s++) all_beta[j*(j+1)*m/2+s]=betaj[s];
  free(Bj_mat);
  free(betaj);
  free(tauj);
  }
 free(betakj);
 free(yj);
 free(muj);
}

/*******************************************************************************************************/


SEXP gibbs_BALVM_naive(SEXP Rtot_draws, SEXP Rknots, SEXP RNknots, SEXP RknotsI, SEXP Rorder,
			SEXP Rn, SEXP Rm, SEXP Rp, SEXP RGrid_length, SEXP Rgrid, SEXP Rdata, SEXP Rab_sigma,
			SEXP Rab_tau, SEXP Reta_start,SEXP Rallbeta_start,SEXP Rdrws_step, SEXP RBurn)
{
/*
* Gibbs sampler for BALVM with naive prior
*/
 int i, j, k, m, p, n, order, grid_length;
 int pn, t, tot_draws, Nknots;
 int drws_step, it_step, eff_draws;
 int save, dens_comp=0;
 int Burn;
 double *data, *ystar, *knots, *knotsI;
 double *eta, *all_beta;
 double *mu, *Ball_mat, *allmu;
 double a_sigma, b_sigma, a_tau, b_tau;
 double *sigma_vec;
 double *grid;
 char f_mu_name[50], f_sig_name[50]="sigma_draws_naive", f_eta_name[50], f_allmu_name[50]="all_mu_naive";

 FILE *out_mu, *out_sigma, *out_eta, *out_allmu, *in_sigma, *in_allmu;

 SEXP Rsigma;

 Nknots=INTEGER(RNknots)[0];
 drws_step=INTEGER(Rdrws_step)[0];
 m=INTEGER(Rm)[0];
 p=INTEGER(Rp)[0];
 order=INTEGER(Rorder)[0];
 n=INTEGER(Rn)[0];
 Burn=INTEGER(RBurn)[0];

 grid_length=INTEGER(RGrid_length)[0];
 tot_draws=INTEGER(Rtot_draws)[0];
 a_sigma=REAL(Rab_sigma)[0];
 b_sigma=REAL(Rab_sigma)[1];
 a_tau=REAL(Rab_tau)[0];
 b_tau=REAL(Rab_tau)[1];
 Ball_mat=(double *) calloc(n*m*p,sizeof(double));
 sigma_vec=(double *) malloc(p*sizeof(double));
 ystar=(double *) malloc(grid_length*p*sizeof(double));
 
 double ***Pijs;
 Pijs=DM3D(n,p,grid_length);

 pn=p*n;
 mu=(double *) malloc(pn*sizeof(double));
 data=(double *) malloc(pn*sizeof(double));
 eta=(double *) malloc(pn*sizeof(double));

  for (i=0; i<pn; i++) {
    data[i]=REAL(Rdata)[i];
    eta[i]=REAL(Reta_start)[i];
    }

 grid=(double *) malloc(grid_length*sizeof(double));
 for (i=0; i<grid_length; i++) grid[i]=REAL(Rgrid)[i];

 all_beta=(double *) malloc(p*(p+1)*m*sizeof(double)/2);
 for (i=0; i<p*(p+1)*m/2; i++) all_beta[i]=REAL(Rallbeta_start)[i];

 knots=(double *) malloc(Nknots*sizeof(double));
 for (i=0; i<Nknots; i++) knots[i]=REAL(Rknots)[i];
 knotsI=(double *) malloc((Nknots-1)*sizeof(double));
 for (i=0; i<(Nknots-1); i++) knotsI[i]=REAL(RknotsI)[i];


 out_sigma=fopen(f_sig_name, "wb");
 out_allmu=fopen(f_allmu_name, "wb");

 GetRNGstate();
 /*starting draws*/
 it_step=0;
 save=0;
for (eff_draws=0; eff_draws<tot_draws; )
 {
    mu_B_fct (n, p, order, m, Nknots, knots, knotsI,  all_beta, eta, Ball_mat, mu, out_allmu, save);
    one_step_gibbs_naive(m, n, p, Ball_mat, data, mu, a_sigma, b_sigma, a_tau, b_tau, all_beta, sigma_vec );
    comp_etaY_grid(n, p, grid_length, order, m, Nknots, knots, knotsI, grid, all_beta, sigma_vec,
		    data, eta, ystar, Pijs, dens_comp);
    it_step=it_step+1;

if (it_step==drws_step)
 {
  eff_draws=eff_draws+1;
  it_step=0;
  if (eff_draws>(Burn-1))
  {
    save=1;
    sprintf(f_mu_name,"%s%05d","mu_nv_draw",eff_draws+1);
    sprintf(f_eta_name,"%s%05d","eta_nv_draw",eff_draws+1);
    out_mu=fopen(f_mu_name, "w");
    out_eta=fopen(f_eta_name, "w");
    for (i=0; i<n; i++)
    {
      for (k=0; k<p; k++)  {fprintf(out_mu, "%f	", mu[k*n+i]); fprintf(out_eta, "%f	", eta[k*n+i]);}
      fprintf(out_mu, "\n"); fprintf(out_eta, "\n");
    }
    fclose(out_mu);
    fclose(out_eta);
    fsafewrite(sigma_vec, sizeof(double),p, out_sigma);
  } else save=0;
 } else save=0;

Rprintf("it_step=%d,  eff_draws=%d\n",it_step, eff_draws);
 }
 PutRNGstate();

 fclose(out_allmu);
 fclose(out_sigma);


 free(data);
 free(grid);
 free(knots);
 free(knotsI);
 free(Ball_mat);
 free(ystar);
 free(mu);
 free(eta);
 free(all_beta);
 Dfree3d(Pijs);

  allmu=(double *) calloc(n*p*(p+1)/2,sizeof(double));
  in_sigma=fopen(f_sig_name, "rb");
  in_allmu=fopen(f_allmu_name, "rb");
  PROTECT(Rsigma = allocMatrix(REALSXP,  p, (tot_draws-Burn)));

 for (t=Burn-1; t<tot_draws-1; t++)
 {
    fsaferead(sigma_vec, sizeof(double), p, in_sigma);
    for(i = 0; i < p; i++) REAL(Rsigma)[(t-Burn+1)*p+i] = sigma_vec[i];
    /*all_mu processing*/
    fsaferead(allmu, sizeof(double), p*(p+1)*n/2, in_allmu);
    sprintf(f_eta_name,"%s%05d","allmu_nv_draw",t+1);
    for (j=0; j<p; j++)
    {
      sprintf(f_mu_name,"%s%s%03d",f_eta_name, "j", j+1);
      out_allmu=fopen(f_mu_name, "w");
      for (i=0; i<n; i++) {
	for (k=0; k<j+1; k++) fprintf(out_allmu, "%f	", allmu[j*(j+1)*n/2+k*n+i]);
	fprintf(out_allmu, "\n");
      }
      fclose(out_allmu);
    }
    /*end all_mu processing*/
   }
  UNPROTECT(1);
  free(sigma_vec);
  free(allmu);
  fclose(in_sigma);
  fclose(in_allmu);
  return (Rsigma);
}

/**************************************************************************/
void getXtWX(double *XtWX, double *X,double *w,int *r,int *c,double *work)
/* forms X'WX as efficiently as possible, where W = diag(w)
   and X is an r by c matrix stored column wise.
   work should be an r-vector (longer is no problem).
*/
{ int i,j;
  double *p,*p1,*p2,*pX0,*pX1,xx;
  pX0=X;
  for (i=0;i< *c;i++) {
    p2 = work + *r;
    for (p=w,p1=work;p1<p2;p++,p1++,pX0++) *p1 = *pX0 * *p;
    for (pX1=X,j=0;j<=i;j++) {
      for (xx=0.0,p=work;p<p2;p++,pX1++) xx += *p * *pX1;
      XtWX[i * *c + j] = XtWX[j * *c + i] = xx;
    }
  }
}
/**************************************************************************/


/*******************************************************************************************************/
/*computing marginal densities estimation */

SEXP jpdf_y1y2(SEXP y1, SEXP y2, SEXP mu11, SEXP sigma1, SEXP mu21, SEXP mu22, SEXP sigma2, SEXP Rn)
{
  int n, nc, nr, i, j, t,s, zr=0;
  double dnsy1, dnsy2;
  dnsy1=0; dnsy2=0;
  SEXP rslt1, rslt2, rslt3, rslt4, ans;

  n=INTEGER(Rn)[0];
  nr=length(y1);
  nc=length(y2);

  PROTECT(rslt1 = allocMatrix(REALSXP, nr, nc));
  PROTECT(rslt2 = allocMatrix(REALSXP, nr, nc));
  PROTECT(rslt3 = allocVector(REALSXP, nr));
  PROTECT(rslt4 = allocVector(REALSXP, nc));
  for (j=0; j<nc; j++){for (i=0; i<nr; i++) REAL(rslt1)[j*nr+i]=0;}


  for (j=0; j<nc; j++)
  {
    for (i=0; i<nr; i++)
    {
      dnsy1=0; dnsy2=0;
      for (s=0; s<n; s++)
      {
	dnsy1=dnsy1+dnorm(REAL(y1)[i], REAL(mu11)[s], REAL(sigma1)[0], zr)/n;
	for (t=0; t<n; t++)
	{
	  dnsy2=dnsy2+dnorm(REAL(y2)[j], REAL(mu21)[s]+REAL(mu22)[t], REAL(sigma2)[0], zr)/(n*n);
	  REAL(rslt1)[j*nr+i]=REAL(rslt1)[j*nr+i]+dnorm(REAL(y1)[i], REAL(mu11)[s], REAL(sigma1)[0], zr)*dnorm(REAL(y2)[j], REAL(mu21)[s]+REAL(mu22)[t], REAL(sigma2)[0], zr)/(n*n);
	}
      }
      REAL(rslt2)[j*nr+i]=REAL(rslt1)[j*nr+i]/(dnsy1*dnsy2);
      REAL(rslt3)[i]=dnsy1;
    }
    REAL(rslt4)[j]=dnsy2;
  }
  PROTECT(ans = allocVector(VECSXP, 4));
  SET_VECTOR_ELT(ans, 0, rslt1);
  SET_VECTOR_ELT(ans, 1, rslt2);
  SET_VECTOR_ELT(ans, 2, rslt3);
  SET_VECTOR_ELT(ans, 3, rslt4);
  UNPROTECT(5);
  return ans;
}

/*******************************************************************************************************/

SEXP jpdf_y2y3(SEXP y2, SEXP y3, SEXP mu21, SEXP mu22, SEXP sigma2, SEXP mu31, SEXP mu32, SEXP mu33, SEXP sigma3, SEXP Rn)
{
  int n, nc, nr, i, j, t,s,r,zr=0;
  double dnsy3, dnsy2;

  SEXP rslt1, rslt2, rslt3, rslt4, ans;

  n=INTEGER(Rn)[0];
  nr=length(y2);
  nc=length(y3);
  dnsy3=0; dnsy2=0;
  PROTECT(rslt1 = allocMatrix(REALSXP, nr, nc));
  PROTECT(rslt2 = allocMatrix(REALSXP, nr, nc));
  PROTECT(rslt3 = allocVector(REALSXP, nr));
  PROTECT(rslt4 = allocVector(REALSXP, nc));
  for (j=0; j<nc; j++){for (i=0; i<nr; i++) REAL(rslt1)[j*nr+i]=0;}


  for (j=0; j<nc; j++)
  {
    for (i=0; i<nr; i++)
    {
      dnsy3=0; dnsy2=0;
      for (s=0; s<n; s++)
      {
	for (t=0; t<n; t++)
	{
	  dnsy2=dnsy2+dnorm(REAL(y2)[i], REAL(mu21)[s]+REAL(mu22)[t], REAL(sigma2)[0], zr)/(n*n);
	  for (r=0; r<n; r++)
	  {
	    dnsy3=dnsy3+dnorm(REAL(y3)[j], REAL(mu31)[s]+REAL(mu32)[t]+REAL(mu33)[r], REAL(sigma3)[0], zr)/(n*n*n);
	    REAL(rslt1)[j*nr+i]=REAL(rslt1)[j*nr+i]+dnorm(REAL(y2)[i], REAL(mu21)[s]+REAL(mu22)[t], REAL(sigma2)[0], zr)*dnorm(REAL(y3)[j], REAL(mu31)[s]+REAL(mu32)[t]+REAL(mu33)[r], REAL(sigma3)[0], zr)/(n*n*n);
	  }
	}
      }
      REAL(rslt2)[j*nr+i]=REAL(rslt1)[j*nr+i]/(dnsy2*dnsy3);
      REAL(rslt3)[i]=dnsy2;
    }
    REAL(rslt4)[j]=dnsy3;
  }
  PROTECT(ans = allocVector(VECSXP, 4));
  SET_VECTOR_ELT(ans, 0, rslt1);
  SET_VECTOR_ELT(ans, 1, rslt2);
  SET_VECTOR_ELT(ans, 2, rslt3);
  SET_VECTOR_ELT(ans, 3, rslt4);
  UNPROTECT(5);
  return ans;
}
/*******************************************************************************************************/

SEXP jpdf_y1y3(SEXP y1, SEXP y3, SEXP mu11, SEXP sigma1, SEXP mu31, SEXP mu32, SEXP mu33, SEXP sigma3, SEXP Rn)
{
  int n, nc, nr, i, j, t,s,r,zr=0;
  double dnsy1, dnsy3;
  dnsy1=0; dnsy3=0;
  SEXP rslt1, rslt2, rslt3, rslt4, ans;
  n=INTEGER(Rn)[0];
  nr=length(y1);
  nc=length(y3);
  PROTECT(rslt1 = allocMatrix(REALSXP, nr, nc));
  PROTECT(rslt2 = allocMatrix(REALSXP, nr, nc));
  PROTECT(rslt3 = allocVector(REALSXP, nr));
  PROTECT(rslt4 = allocVector(REALSXP, nc));
  for (j=0; j<nc; j++){for (i=0; i<nr; i++) REAL(rslt1)[j*nr+i]=0;}


  for (j=0; j<nc; j++)
  {
    for (i=0; i<nr; i++)
    {
      dnsy1=0; dnsy3=0;
      for (s=0; s<n; s++)
      {
	dnsy1=dnsy1+dnorm(REAL(y1)[i], REAL(mu11)[s], REAL(sigma1)[0], zr)/n;
	for (t=0; t<n; t++)
	{
	  for (r=0; r<n; r++)
	  {
	    dnsy3=dnsy3+dnorm(REAL(y3)[j], REAL(mu31)[s]+REAL(mu32)[t]+REAL(mu33)[r], REAL(sigma3)[0], zr)/(n*n*n);
	    REAL(rslt1)[j*nr+i]=REAL(rslt1)[j*nr+i]+dnorm(REAL(y1)[i], REAL(mu11)[s], REAL(sigma1)[0], zr)*dnorm(REAL(y3)[j], REAL(mu31)[s]+REAL(mu32)[t]+REAL(mu33)[r], REAL(sigma3)[0], zr)/(n*n*n);
	  }
	}
      }
      REAL(rslt2)[j*nr+i]=REAL(rslt1)[j*nr+i]/(dnsy1*dnsy3);
      REAL(rslt3)[i]=dnsy1;
    }
    REAL(rslt4)[j]=dnsy3;
  }
  PROTECT(ans = allocVector(VECSXP, 4));
  SET_VECTOR_ELT(ans, 0, rslt1);
  SET_VECTOR_ELT(ans, 1, rslt2);
  SET_VECTOR_ELT(ans, 2, rslt3);
  SET_VECTOR_ELT(ans, 3, rslt4);
  UNPROTECT(5);
  return ans;
}

/*******************************************************************************************************/

SEXP cpdf_y2y3givy1(SEXP y1, SEXP y2, SEXP y3, SEXP mu11, SEXP sigma1, SEXP mu21, SEXP mu22, SEXP sigma2, SEXP mu31, SEXP mu32, SEXP mu33, SEXP sigma3, SEXP Rn)
{
  int n, nc, nr, i, j, t,s,r,zr=0;
  double dnsy1;
  SEXP rslt;
  n=INTEGER(Rn)[0];
  nr=length(y2);
  nc=length(y3);
  PROTECT(rslt = allocMatrix(REALSXP, nr, nc));
  for (j=0; j<nc; j++){for (i=0; i<nr; i++) REAL(rslt)[j*nr+i]=0;}
  dnsy1=0;
  for (s=0; s<n; s++) dnsy1=dnsy1+dnorm(REAL(y1)[0],REAL(mu11)[s], REAL(sigma1)[0], zr);

  for (j=0; j<nc; j++)
  {
    for (i=0; i<nr; i++)
    {
      for (s=0; s<n; s++)
      {
	for (t=0; t<n; t++)
	{
	  for (r=0; r<n; r++) REAL(rslt)[j*nr+i]=REAL(rslt)[j*nr+i]+dnorm(REAL(y1)[0],REAL(mu11)[s], REAL(sigma1)[0], zr)*
	    dnorm(REAL(y2)[i],REAL(mu21)[s]+REAL(mu22)[t], REAL(sigma2)[0], zr)*dnorm(REAL(y3)[j], REAL(mu31)[s]+REAL(mu32)[t]+REAL(mu33)[r], REAL(sigma3)[0], zr)/(dnsy1*n*n);
	}
      }
    }
  }
  UNPROTECT(1);
  return rslt;
}

/*******************************************************************************************************/

SEXP wrong_cpdf_y2y3givy1(SEXP y1, SEXP y2, SEXP y3, SEXP mu1, SEXP sigma1, SEXP mu2, SEXP sigma2, SEXP mu3, SEXP sigma3, SEXP Rn)
{
  int n, nc, nr, i, j,s,zr=0;
  double dnsy1;
  SEXP rslt;
  n=INTEGER(Rn)[0];
  nr=length(y2);
  nc=length(y3);
  PROTECT(rslt = allocMatrix(REALSXP, nr, nc));
  for (j=0; j<nc; j++){for (i=0; i<nr; i++) REAL(rslt)[j*nr+i]=0;}
  dnsy1=0;
  for (s=0; s<n; s++) dnsy1=dnsy1+dnorm(REAL(y1)[0],REAL(mu1)[s], REAL(sigma1)[0], zr);

  for (j=0; j<nc; j++)
  {
    for (i=0; i<nr; i++)
    {
      for (s=0; s<n; s++)
	REAL(rslt)[j*nr+i]=REAL(rslt)[j*nr+i]+dnorm(REAL(y1)[0],REAL(mu1)[s], REAL(sigma1)[0], zr)*
	    dnorm(REAL(y2)[i],REAL(mu2)[s], REAL(sigma2)[0], zr)*dnorm(REAL(y3)[j], REAL(mu3)[s], REAL(sigma3)[0], zr)/(dnsy1);

      }
   }
  UNPROTECT(1);
  return rslt;
}

/*******************************************************************************************************/

SEXP wrong_jpdf_y2y3(SEXP y2, SEXP y3, SEXP mu2, SEXP sigma2, SEXP mu3, SEXP sigma3, SEXP Rn)
{
  int n, nc, nr, i, j, s,zr=0;
  double dnsy3, dnsy2;

  SEXP rslt1;

  n=INTEGER(Rn)[0];
  nr=length(y2);
  nc=length(y3);
  dnsy3=0; dnsy2=0;
  PROTECT(rslt1 = allocMatrix(REALSXP, nr, nc));
  for (j=0; j<nc; j++){for (i=0; i<nr; i++) REAL(rslt1)[j*nr+i]=0;}

  for (j=0; j<nc; j++)
  {
    for (i=0; i<nr; i++)
    {
      dnsy3=0; dnsy2=0;
      for (s=0; s<n; s++)
      {
	    dnsy3=dnsy3+dnorm(REAL(y3)[j], REAL(mu3)[s], REAL(sigma3)[0], zr)/(n);
	    dnsy2=dnsy2+dnorm(REAL(y2)[i], REAL(mu2)[s], REAL(sigma2)[0], zr)/(n);
	    REAL(rslt1)[j*nr+i]=REAL(rslt1)[j*nr+i]+dnorm(REAL(y2)[i], REAL(mu2)[s], REAL(sigma2)[0], zr)*dnorm(REAL(y3)[j], REAL(mu3)[s], REAL(sigma3)[0], zr)/n;
      }
    }
  }
  UNPROTECT(1);
  return rslt1;
}

/*******************************************************************************************************/

SEXP wrong_jpdf_y1y3(SEXP y1, SEXP y3, SEXP mu1, SEXP sigma1, SEXP mu3, SEXP sigma3, SEXP Rn)
{
  int n, nc, nr, i, j, s,zr=0;
  double dnsy3, dnsy1;

  SEXP rslt1;

  n=INTEGER(Rn)[0];
  nr=length(y1);
  nc=length(y3);
  dnsy3=0; dnsy1=0;
  PROTECT(rslt1 = allocMatrix(REALSXP, nr, nc));
  for (j=0; j<nc; j++){for (i=0; i<nr; i++) REAL(rslt1)[j*nr+i]=0;}


  for (j=0; j<nc; j++)
  {
    for (i=0; i<nr; i++)
    {
      dnsy3=0; dnsy1=0;
      for (s=0; s<n; s++)
      {
	    dnsy3=dnsy3+dnorm(REAL(y3)[j], REAL(mu3)[s], REAL(sigma3)[0], zr)/(n);
	    dnsy1=dnsy1+dnorm(REAL(y1)[i], REAL(mu1)[s], REAL(sigma1)[0], zr)/(n);
	    REAL(rslt1)[j*nr+i]=REAL(rslt1)[j*nr+i]+dnorm(REAL(y1)[i], REAL(mu1)[s], REAL(sigma1)[0], zr)*dnorm(REAL(y3)[j], REAL(mu3)[s], REAL(sigma3)[0], zr)/n;
      }
    }

  }
  UNPROTECT(1);
  return rslt1;
}

/*******************************************************************************************************/

SEXP wrong_jpdf_y1y2(SEXP y1, SEXP y2, SEXP mu1, SEXP sigma1, SEXP mu2, SEXP sigma2, SEXP Rn)
{
  int n, nc, nr, i, j, s,zr=0;
  double dnsy2, dnsy1;

  SEXP rslt1;

  n=INTEGER(Rn)[0];
  nr=length(y1);
  nc=length(y2);
  dnsy2=0; dnsy1=0;
  PROTECT(rslt1 = allocMatrix(REALSXP, nr, nc));
  for (j=0; j<nc; j++){for (i=0; i<nr; i++) REAL(rslt1)[j*nr+i]=0;}


  for (j=0; j<nc; j++)
  {
    for (i=0; i<nr; i++)
    {
      dnsy2=0; dnsy1=0;
      for (s=0; s<n; s++)
      {
	    dnsy2=dnsy2+dnorm(REAL(y2)[j], REAL(mu2)[s], REAL(sigma2)[0], zr)/(n);
	    dnsy1=dnsy1+dnorm(REAL(y1)[i], REAL(mu1)[s], REAL(sigma1)[0], zr)/(n);
	    REAL(rslt1)[j*nr+i]=REAL(rslt1)[j*nr+i]+dnorm(REAL(y1)[i], REAL(mu1)[s], REAL(sigma1)[0], zr)*dnorm(REAL(y2)[j], REAL(mu2)[s], REAL(sigma2)[0], zr)/n;
      }
    }
  }
  UNPROTECT(1);
  return rslt1;
}


/* void wjpdf_y1y2(double* y1, int ly1, double* y2, int ly2, double* mu1, double sigma1, double* mu2, double sigma2, int n, double **dens_mat)
{
  int i, j, s,zr=0;
  double dnsy2, dnsy1;
  dnsy2=0; dnsy1=0;

  for (j=0; j<ly2; j++)
  {
    for (i=0; i<ly1; i++)
    {
      dnsy2=0; dnsy1=0;
      for (s=0; s<n; s++)
      {
	    dnsy2=dnsy2+dnorm(y2[j], mu2[s], sigma2, zr)/(n);
	    dnsy1=dnsy1+dnorm(y1[i], mu1[s], sigma1[0], zr)/(n);
	    dens_mat[i][j]=dens_mat[i][j]+dnorm(y1[i], mu1[s], sigma1, zr)*dnorm(y2[j], mu2[s], sigma2, zr)/n;
      }
    }
  }
}



SEXP mean_density(SEXP Ry1, SEXP Ry2, SEXP Ry3, SEXP mu1, SEXP sigma1, SEXP mu2, SEXP sigma2, SEXP mu3, SEXP sigma3, SEXP Rn, SEXP Rburn, SEXP Rtot_draws)
{
 int ly1, ly2, ly3;
 ly1=length(y1);
 ly2=length(y2);
 ly3=length(y3);
 int n, i;
 n=INTEGER(Rn)[0];
 double ***Pijs;
 Pijs=DM3D(ly1,ly2,INTEGER(Rtot_draws)[0]-INTEGER(Rburn)[0]);
 double **dens_mat;
 dens_mat=DM2D(ly1,ly2);
for (i=0; i<INTEGER(Rtot_draws)[0]-INTEGER(Rburn)[0]; i++)

 Dfree3d(Pijs);
}*/

/*******************************************************************************************************/

SEXP L1_KL_dist_hard(SEXP Ry1, SEXP Ry2, SEXP Ry3, SEXP Rmu1, SEXP Rsigma1, SEXP Rmu2, SEXP Rsigma2, SEXP Rmu3, SEXP Rsigma3,
		 SEXP Tmu1, SEXP Tsigma1, SEXP Tmu2, SEXP Tsigma2, SEXP Tmu3, SEXP Tsigma3, SEXP Rn)
{
  int n,s,zr=0;
  int ly1, ly2, ly3, i,j,k;
  double dnsym, dnsyt;
  SEXP ans;
  ly1=length(Ry1);ly2=length(Ry2);ly3=length(Ry3);

  n=INTEGER(Rn)[0];
  dnsym=0; dnsyt=0;
  PROTECT(ans = allocVector(REALSXP, 2));
  REAL(ans)[0]=0;
  REAL(ans)[1]=0;

  if ((ly1==ly2)&(ly2==ly3))
  {
    for (i=0; i<ly1; i++)
    {
     dnsym=0; dnsyt=0;
     for (s=0; s<n; s++)
      {
	    dnsym=dnsym+dnorm(REAL(Ry1)[i], REAL(Rmu1)[s], REAL(Rsigma1)[0], zr)*
			dnorm(REAL(Ry2)[i], REAL(Rmu2)[s], REAL(Rsigma2)[0], zr)*
			dnorm(REAL(Ry3)[i], REAL(Rmu3)[s], REAL(Rsigma3)[0], zr)/n;
      }
    double mu1, mu2, mu3, s1,s2, seps;
    mu1=55; mu2=10; mu3=80; s1=5; s2=2.25; seps=1;
      dnsyt=dnorm(REAL(Ry1)[i], mu1, sqrt(pow(s1,2)+pow(seps, 2)), zr)*
    		dnorm(REAL(Ry2)[i], mu2, sqrt(pow(s2,2)+pow(seps, 2)), zr)*
    		dnorm(REAL(Ry3)[i], mu3+2*REAL(Ry2)[i]*exp(REAL(Ry1)[i]-50)/(1+exp(REAL(Ry1)[i]-50)), seps, zr);
     REAL(ans)[0]=REAL(ans)[0]+fabs(dnsym-dnsyt)/(ly1);
     Rprintf("dnsym=%f,     dnsyt=%f\n", dnsym, dnsyt);
     REAL(ans)[1]=REAL(ans)[1]+0.5*dnsym*log(dnsym/dnsyt)/(ly1)+0.5*dnsyt*log(dnsyt/dnsym)/(ly1);
    }
  } else
  {
    for (i=0; i<ly1; i++)
    {
      for (j=0; j<ly2; j++)
      {
	for (k=0; k<ly3; k++)
	{
	  dnsym=0; dnsym=0;
	  for (s=0; s<n; s++)
	  {
		dnsym=dnsym+dnorm(REAL(Ry1)[i], REAL(Rmu1)[s], REAL(Rsigma1)[0], zr)*
			    dnorm(REAL(Ry2)[j], REAL(Rmu2)[s], REAL(Rsigma2)[0], zr)*
			    dnorm(REAL(Ry3)[k], REAL(Rmu3)[s], REAL(Rsigma3)[0], zr)/n;
	  }
    double mu1, mu2, mu3, s1,s2, seps;
    mu1=55; mu2=10; mu3=80; s1=5; s2=2.25; seps=1;
      dnsyt=dnorm(REAL(Ry1)[i], mu1, sqrt(pow(s1,2)+pow(seps, 2)), zr)*
    		dnorm(REAL(Ry2)[j], mu2, sqrt(pow(s2,2)+pow(seps, 2)), zr)*
    		dnorm(REAL(Ry3)[k], mu3+2*REAL(Ry2)[i]*exp(REAL(Ry1)[i]-50)/(1+exp(REAL(Ry1)[i]-50)), seps, zr);
	  REAL(ans)[0]=REAL(ans)[0]+fabs(dnsym-dnsyt)/(ly1*ly2*ly3);
	  REAL(ans)[1]=REAL(ans)[1]+0.5*dnsym*log(dnsym/dnsyt)/(ly1*ly2*ly3)+0.5*dnsyt*log(dnsyt/dnsym)/(ly1*ly2*ly3);
	}
      }
    }
  }
  UNPROTECT(1);
  return ans;
}

/*******************************************************************************************************/

SEXP L1_KL_dist_simple(SEXP Ry1, SEXP Ry2, SEXP Ry3, SEXP Rmu1, SEXP Rsigma1, SEXP Rmu2, SEXP Rsigma2, SEXP Rmu3, SEXP Rsigma3, SEXP Rn)
{
  int n,s,zr=0;
  int ly1, ly2, ly3, i,j,k;
  double dnsym, dnsyt;
  double w11=0.6,w12,w21=0.9,w22,c=2.5, s1=1,s2=0.8,s3=0.7,seps=1, mu1=-1,mu2=4,mu3=3;
  SEXP ans;
  w12=1-w11; w22=1-w21;
  ly1=length(Ry1);ly2=length(Ry2);ly3=length(Ry3);

  n=INTEGER(Rn)[0];
  dnsym=0; dnsyt=0;
  PROTECT(ans = allocVector(REALSXP, 2));
  REAL(ans)[0]=0;
  REAL(ans)[1]=0;

  if ((ly1==ly2)&(ly2==ly3))
  {
    for (i=0; i<ly1; i++)
    {
     dnsym=0; dnsyt=0;
     for (s=0; s<n; s++)
      {
	    dnsym=dnsym+dnorm(REAL(Ry1)[i], REAL(Rmu1)[s], REAL(Rsigma1)[0], zr)*
			dnorm(REAL(Ry2)[i], REAL(Rmu2)[s], REAL(Rsigma2)[0], zr)*
			dnorm(REAL(Ry3)[i], REAL(Rmu3)[s], REAL(Rsigma3)[0], zr)/n;
      }
      dnsyt=(w11*dnorm(REAL(Ry1)[i], mu1, s1, zr)+w12*dnorm(REAL(Ry1)[i], mu2, s2, zr))*
	    (w21*dnorm(REAL(Ry2)[i], mu1, s1, zr)+w22*dnorm(REAL(Ry2)[i], mu3, s3, zr))*
	    dnorm(REAL(Ry3)[i], REAL(Ry1)[i]+c*REAL(Ry2)[i], seps, zr);
      REAL(ans)[0]=REAL(ans)[0]+fabs(dnsym-dnsyt)/(ly1);
     /*Rprintf("dnsym=%f,     dnsyt=%f\n", dnsym, dnsyt);*/
     REAL(ans)[1]=REAL(ans)[1]+0.5*dnsym*log(dnsym/dnsyt)/(ly1)+0.5*dnsyt*log(dnsyt/dnsym)/(ly1);
    }
  } else
  {
    for (i=0; i<ly1; i++)
    {
      for (j=0; j<ly2; j++)
      {
	for (k=0; k<ly3; k++)
	{
	  dnsym=0; dnsym=0;
	  for (s=0; s<n; s++)
	  {
		dnsym=dnsym+dnorm(REAL(Ry1)[i], REAL(Rmu1)[s], REAL(Rsigma1)[0], zr)*
			    dnorm(REAL(Ry2)[j], REAL(Rmu2)[s], REAL(Rsigma2)[0], zr)*
			    dnorm(REAL(Ry3)[k], REAL(Rmu3)[s], REAL(Rsigma3)[0], zr)/n;
	  }
    dnsyt=(w11*dnorm(REAL(Ry1)[i], mu1, sqrt(pow(s1,2)+pow(seps, 2)), zr)+w12*dnorm(REAL(Ry1)[i], mu2, sqrt(pow(s2,2)+pow(seps, 2)), zr))*
  		(w21*dnorm(REAL(Ry2)[i], mu1, sqrt(pow(s1,2)+pow(seps, 2)), zr)+w22*dnorm(REAL(Ry2)[i], mu3, sqrt(pow(s3,2)+pow(seps, 2)), zr))*
  		dnorm(REAL(Ry3)[i], REAL(Ry1)[i]+c*REAL(Ry2)[i], seps, zr);
	  REAL(ans)[0]=REAL(ans)[0]+fabs(dnsym-dnsyt)/(ly1*ly2*ly3);
	  REAL(ans)[1]=REAL(ans)[1]+0.5*dnsym*log(dnsym/dnsyt)/(ly1*ly2*ly3)+0.5*dnsyt*log(dnsyt/dnsym)/(ly1*ly2*ly3);
	  Rprintf("dnsym=%f,     dnsyt=%f\n", dnsym, dnsyt);
	}
      }
    }
  }
  UNPROTECT(1);
  return ans;
}


/*******************************************************************************************************/

SEXP mut_infos(SEXP Rv1, SEXP Rv2, SEXP Rv3, SEXP Rev1, SEXP Rsigma1, SEXP Rev2, SEXP Rsigma2, SEXP Rev3, SEXP Rsigma3, SEXP Rn)
{
  int n, i, j, k, s,zr=0;
  double dnsv1, dnsv2, dnsv3, v1v2dns, v2v3dns, v1v3dns;
  int lv1,lv2,lv3;
  SEXP ans;

  n=INTEGER(Rn)[0];

  lv1=length(Rv1);
  lv2=length(Rv2);
  lv3=length(Rv3);

  dnsv3=0; dnsv2=0; dnsv1=0;
  v1v2dns=0; v2v3dns=0; v1v3dns=0;

  PROTECT(ans = allocVector(REALSXP, 3));
  for (j=0; j<3; j++) REAL(ans)[j]=0;


  if ((lv1==lv2)&(lv2==lv3))
  {
    for (i=0; i<lv1; i++)
    {
      dnsv3=0; dnsv2=0; dnsv1=0;
      for (s=0; s<n; s++)
      {
	    dnsv3=dnsv3+dnorm(REAL(Rv3)[i], REAL(Rev3)[s], REAL(Rsigma3)[0], zr)/(n);
	    dnsv2=dnsv2+dnorm(REAL(Rv2)[i], REAL(Rev2)[s], REAL(Rsigma2)[0], zr)/(n);
	    dnsv1=dnsv1+dnorm(REAL(Rv1)[i], REAL(Rev1)[s], REAL(Rsigma1)[0], zr)/(n);
	    v1v2dns= v1v2dns+dnorm(REAL(Rv1)[i], REAL(Rev1)[s], REAL(Rsigma1)[0], zr)*dnorm(REAL(Rv2)[i], REAL(Rev2)[s], REAL(Rsigma2)[0], zr)/n;
	    v2v3dns= v2v3dns+dnorm(REAL(Rv2)[i], REAL(Rev2)[s], REAL(Rsigma2)[0], zr)*dnorm(REAL(Rv3)[i], REAL(Rev3)[s], REAL(Rsigma3)[0], zr)/n;
	    v1v3dns= v1v3dns+dnorm(REAL(Rv1)[i], REAL(Rev1)[s], REAL(Rsigma1)[0], zr)*dnorm(REAL(Rv3)[i], REAL(Rev3)[s], REAL(Rsigma3)[0], zr)/n;
      }
      REAL(ans)[0]=REAL(ans)[0]+v1v2dns*log(v1v2dns/(dnsv1*dnsv2))/lv1;
      REAL(ans)[1]=REAL(ans)[1]+v2v3dns*log(v2v3dns/(dnsv2*dnsv3))/lv1;
      REAL(ans)[2]=REAL(ans)[2]+v1v3dns*log(v1v3dns/(dnsv1*dnsv3))/lv1;
  }
  } else {
  for (i=0; i<lv1; i++)
  {
    for (j=0; j<lv2; j++)
    {
       for (k=0; k<lv3; k++)
       {

      dnsv3=0; dnsv2=0; dnsv1=0;
      for (s=0; s<n; s++)
      {
	    dnsv3=dnsv3+dnorm(REAL(Rv3)[k], REAL(Rev3)[s], REAL(Rsigma3)[0], zr)/(n);
	    dnsv2=dnsv2+dnorm(REAL(Rv2)[j], REAL(Rev2)[s], REAL(Rsigma2)[0], zr)/(n);
	    dnsv1=dnsv1+dnorm(REAL(Rv1)[i], REAL(Rev1)[s], REAL(Rsigma1)[0], zr)/(n);
	    v1v2dns= v1v2dns+dnorm(REAL(Rv1)[i], REAL(Rev1)[s], REAL(Rsigma1)[0], zr)*dnorm(REAL(Rv2)[j], REAL(Rev2)[s], REAL(Rsigma2)[0], zr)/n;
	    v2v3dns= v2v3dns+dnorm(REAL(Rv2)[j], REAL(Rev2)[s], REAL(Rsigma2)[0], zr)*dnorm(REAL(Rv3)[k], REAL(Rev3)[s], REAL(Rsigma3)[0], zr)/n;
	    v1v3dns= v1v3dns+dnorm(REAL(Rv1)[i], REAL(Rev1)[s], REAL(Rsigma1)[0], zr)*dnorm(REAL(Rv3)[k], REAL(Rev3)[s], REAL(Rsigma3)[0], zr)/n;
      }
      REAL(ans)[0]=REAL(ans)[0]+v1v2dns*log(v1v2dns/(dnsv1*dnsv2))/lv1*lv2*lv3;
      REAL(ans)[1]=REAL(ans)[1]+v2v3dns*log(v2v3dns/(dnsv2*dnsv3))/lv1*lv2*lv3;
      REAL(ans)[2]=REAL(ans)[2]+v1v3dns*log(v1v3dns/(dnsv1*dnsv3))/lv1*lv2*lv3;
    }
  }
  }
  }
  UNPROTECT(1);
  return ans;
}

/*******************************************************************************************************/

SEXP triv_DPMNdens(SEXP Rv1, SEXP Rv2, SEXP Rv3, SEXP Rmu1, SEXP Rmu2, SEXP Rmu3, SEXP Rsigvec, SEXP Rsimple)
{
  int i,s;
  int l1;
  int zr=0;
  int lv1;
  double v1v2v3dns, dnsyt;
  SEXP ans;
  l1=length(Rmu1);  /*=length(Rmu2)=length(Rmu3);*/
  lv1=length(Rv1); /*=length(Rv2)=length(Rv3);*/
  v1v2v3dns=0; dnsyt=0;

  PROTECT(ans = allocVector(REALSXP, 2));
  REAL(ans)[0]=0;
  REAL(ans)[1]=0;

  for (i=0; i<lv1; i++)
   {
     v1v2v3dns=0;
	   for (s=0; s<l1; s++)
	   {
  	  v1v2v3dns=v1v2v3dns+dnorm(REAL(Rv1)[i], REAL(Rmu1)[s], REAL(Rsigvec)[0], zr)*
			dnorm(REAL(Rv2)[i], REAL(Rmu2)[s], REAL(Rsigvec)[1], zr)*
			dnorm(REAL(Rv3)[i], REAL(Rmu3)[s], REAL(Rsigvec)[2], zr)/l1;
			}
    if (INTEGER(Rsimple)[0]==1)
    {
      double w11=0.6,w12,w21=0.9,w22,c=2.5, s1=1,s2=0.8,s3=0.7,seps=1, mu1=-1,mu2=4,mu3=3;
      w12=1-w11; w22=1-w21;
      dnsyt=(w11*dnorm(REAL(Rv1)[i], mu1, sqrt(pow(s1,2)+pow(seps, 2)), zr)+w12*dnorm(REAL(Rv1)[i], mu2, sqrt(pow(s2,2)+pow(seps, 2)), zr))*
    		(w21*dnorm(REAL(Rv2)[i], mu1, sqrt(pow(s1,2)+pow(seps, 2)), zr)+w22*dnorm(REAL(Rv2)[i], mu3, sqrt(pow(s3,2)+pow(seps, 2)), zr))*
    		dnorm(REAL(Rv3)[i], REAL(Rv1)[i]+c*REAL(Rv2)[i], seps, zr);
	  } else
    {
    double mu1, mu2, mu3, s1,s2, seps;
    mu1=55; mu2=10; mu3=80; s1=5; s2=2.25; seps=1;
      dnsyt=dnorm(REAL(Rv1)[i], mu1, sqrt(pow(s1,2)+pow(seps, 2)), zr)*
    		dnorm(REAL(Rv2)[i], mu2, sqrt(pow(s2,2)+pow(seps, 2)), zr)*
    		dnorm(REAL(Rv3)[i], mu3+2*REAL(Rv2)[i]*exp(REAL(Rv1)[i]-50)/(1+exp(REAL(Rv1)[i]-50)), seps, zr);
    }
  		REAL(ans)[0]=REAL(ans)[0]+fabs(v1v2v3dns-dnsyt)/(lv1);
   	  if (((v1v2v3dns/dnsyt)>1e-20)&((dnsyt/v1v2v3dns)>1e-20)) REAL(ans)[1]=REAL(ans)[1]+0.5*v1v2v3dns*log(v1v2v3dns/dnsyt)/(lv1)+0.5*dnsyt*log(dnsyt/v1v2v3dns)/(lv1);
		}
 UNPROTECT(1);
 return ans;
}


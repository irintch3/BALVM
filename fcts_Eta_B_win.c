/*  Routines for manipulating B-splines.  These are intended for use with
 *  S or S-PLUS or R.
 *
 *     Copyright (C) 1998 Douglas M. Bates and William N. Venables.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 *
 * These functions are distributed in the hope that they will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the
 * GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 *
 * The routines are loosely based on the pseudo-code in Schumacher (Wiley,
 * 1981) and the CMLIB library DBSPLINES.
 */

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

typedef struct {
    int order,			/* order of the spline */
	ordm1,			/* order - 1 (3 for cubic splines) */
	nknots,			/* number of knots */
	curs,			/* current position in knots vector */
	boundary;		/* must have knots[curs] <= x < knots[curs+1] */
				/* except for the boundary case */

    double *ldel,		/* differences from knots on the left */
	*rdel,			/* differences from knots on the right */
	*knots,			/* knot vector */
	*coeff,			/* coefficients */
	*a;			/* scratch array */
} splPTR;


/* set sp->curs to the index of the first knot position > x.
   Special handling for x == sp->knots[sp->nknots - sp-order + 1] */

/*************************************************************************************************************************/
static int
set_cursor(splPTR *sp, double x)
{
    int i;
    /* don't assume x's are sorted */

    sp->curs = -1; /* Wall */
    sp->boundary = 0;
    for (i = 0; i < sp->nknots; i++) {
	if (sp->knots[i] >= x) sp->curs = i;
	if (sp->knots[i] > x) break;
    }
    if (sp->curs > sp->nknots - sp->order) {
	int lastLegit = sp->nknots - sp->order;
	if (x == sp->knots[lastLegit]) {
	    sp->boundary = 1; sp->curs = lastLegit;
	}
    }
    return sp->curs;
}

/*************************************************************************************************************************/
static void
diff_table(splPTR *sp, double x, int ndiff)
{
  int i;
  for (i = 0; i < ndiff; i++) {
      sp->rdel[i] = sp->knots[sp->curs + i] - x;
      sp->ldel[i] = x - sp->knots[sp->curs - (i + 1)];
  }
}

/******************************************************************************/


void print_matrix( char* desc, int m, int n, double* a, int lda )
{
    int i, j;
    Rprintf( " %s\n", desc );
    for( i = 0; i < m; i++ ) {
            for( j = 0; j < n; j++ ) Rprintf( " %6.9f", a[i+j*lda] );
            Rprintf( "\n" );
    }
}

/*************************************************************************************************************************/
/* fast evaluation of basis functions */
static void
basis_funcs(splPTR *sp, double x, double *b)
{
    int j, r;
    double saved, term;

    diff_table(sp, x, sp->ordm1);
    b[0] = 1.;
    for (j = 1; j <= sp->ordm1; j++) {
	saved = 0.;
	for (r = 0; r < j; r++) {
	    term = b[r]/(sp->rdel[r] + sp->ldel[j - 1 - r]);
	    b[r] = saved + sp->rdel[r] * term;
	    saved = sp->ldel[j - 1 - r] * term;
	}
	b[j] = saved;
    }
}

/*************************************************************************************************************************/
/* "slow" evaluation of (derivative of) basis functions */
static double
evaluate(splPTR *sp, double x, int nder)
{
    register double *lpt, *rpt, *apt, *ti = sp->knots + sp->curs;
    int inner, outer = sp->ordm1;

    if (sp->boundary && nder == sp->ordm1) { /* value is arbitrary */
	return 0.0;
    }
    while(nder--) {
	for(inner = outer, apt = sp->a, lpt = ti - outer; inner--; apt++, lpt++)
	    *apt = outer * (*(apt + 1) - *apt)/(*(lpt + outer) - *lpt);
	outer--;
    }
    diff_table(sp, x, outer);
    while(outer--)
	for(apt = sp->a, lpt = sp->ldel + outer, rpt = sp->rdel, inner = outer + 1;
	    inner--; lpt--, rpt++, apt++)
	    *apt = (*(apt + 1) * *lpt + *apt * *rpt)/(*rpt + *lpt);
    return sp->a[0];
}


/*************************************************************************************************************************/

void spline_basis(double *knots, int nk, int order, double *xvals, int nx, int *derivs, int nd, int *offsets, double *val)
{
/* evaluate the non-zero B-spline basis functions (or their derivatives)
 * at xvals.  */
    int i, j;
    splPTR *sp;
    sp = (splPTR*) calloc(1, sizeof(splPTR));

    /*kk = REAL(knots);*/
    /*xx = REAL(xvals);*/
    /*ders = INTEGER(derivs);*/

    /* fill sp : */
    sp->order = order;
    sp->ordm1 = sp->order - 1;
    sp->rdel = (double *) calloc(sp->ordm1, sizeof(double));
    sp->ldel = (double *) calloc(sp->ordm1, sizeof(double));
    sp->knots = knots; sp->nknots = nk;
    sp->a = (double *) calloc(sp->order, sizeof(double));
    

    for(i = 0; i < nx; i++) {
	set_cursor(sp, xvals[i]);
	offsets[i] = j = sp->curs - sp->order;
	if (j < 0 || j > nk) {
	    for (j = 0; j < sp->order; j++) {
		val[i * sp->order + j] = HIGH;
	    }
	} else {
	    if (derivs[i % nd] > 0) { /* slow method for derivatives */
		int ii;
		for(ii = 0; ii < sp->order; ii++) {
		    for(j = 0; j < sp->order; j++) sp->a[j] = 0;
		    sp->a[ii] = 1;
		    val[i * sp->order + ii] =
			evaluate(sp, xvals[i], derivs[i % nd]);
		}
	    } else {		/* fast method for value */
		basis_funcs(sp, xvals[i], val + i * sp->order);
	    }
	}
    }
free(sp->rdel);
free(sp->ldel);
free(sp->a);
free(sp);
}
/*******************************************************************/
void subvec(double *vec1, int a, int  b, double *resvec)
{
  int i=0, length2=b-a;
  for(i = 0; i < length2; i++) resvec[i]=vec1[a+i];
}

/*******************************************************************/

void all_mu_comp(double *eta, int p, int order, int m, int nknots, double *knots, double *knotsI, double *all_beta, double *all_mu) 
/*all_mu should be an array of zeros*/
{
/* 
m is the degree of freedom of splines and number of components in beta_{kj}
knots and knots is for splines without intercept (but the first coordinates is the intercept function and should be excluded)
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

/************************************************************************/
/* A  version using Walker's alias method, based on Alg 3.13B in
   Ripley (1987).
 */
#define SMALL 10000
static void
walker_ProbSampleReplace(int n, double *p, int *a, int nans, int *ans)
{
    double *q, rU;
    int i, j, k;
    int *HL, *H, *L;


    /* Create the alias tables.
       The idea is that for HL[0] ... L-1 label the entries with q < 1
       and L ... H[n-1] label those >= 1.
       By rounding error we could have q[i] < 1. or > 1. for all entries.
     */

	/* Slow enough anyway not to risk overflow */
	HL = Calloc(n, int);
	q = Calloc(n, double);

    H = HL - 1; L = HL + n;
    for (i = 0; i < n; i++) {
	q[i] = p[i] * n;
	if (q[i] < 1.) *++H = i; else *--L = i;
    }
    if (H >= HL && L < HL + n) { /* So some q[i] are >= 1 and some < 1 */
	for (k = 0; k < n - 1; k++) {
	    i = HL[k];
	    j = *L;
	    a[i] = j;
	    q[j] += q[i] - 1;
	    if (q[j] < 1.) L++;
	    if(L >= HL + n) break; /* now all are >= 1 */
	}
    }
    for (i = 0; i < n; i++) q[i] += i;

    /* generate sample */
    for (i = 0; i < nans; i++) {
	rU = unif_rand() * n;
	k = (int) rU;
	ans[i] = (rU < q[k]) ? k+1 : a[k]+1;
    }
    if(n > SMALL) {
	Free(HL);
	Free(q);
    }
}

/***********************************************************************
dyn.load("fcts_Eta_B_win.so")
Rtot_draws <- as.integer(18)#0000
Rdens=0
Rdrws.step <- 1

mu.all.my <- .Call("gibbs_BALVM_naive", Rtot_draws, Aknots, as.integer(nknots), AknotsI, as.integer(ord),
                   as.integer(n), as.integer(m), as.integer(p), as.integer(S), grid.eta, data, c(1,0.05),
                   c(1,0.05), eta, Rallbeta, as.integer(Rdens),as.integer(Rdrws.step) )

dyn.unload("fcts_Eta_B_win.so")*/

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

/* free the data */
free(**pa);
/* free the row pointers */
free(*pa);
/* free the 3D pointers */
free(pa);
}
/******************************************************************************/
double **DM2D(int n,int row) {
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

free(*pa);
free(pa);
}

/********************************/

long fsafewrite(ptr,size,n,stream)
     double *ptr;size_t size;long n;FILE *stream;

{ long i,j,k=0L;
  for (i=0;i<(n/32000L);i++)
  k+=fwrite(ptr+i*32000L,size,(size_t)32000L,stream);
  j=n%32000L;
  k+=fwrite(ptr+i*32000L,size,(size_t)j,stream);
  return(k);
}

long fsaferead(ptr,size,n,stream)
     double *ptr;size_t size;long n;FILE *stream;

{ long i,j=32000L,k=0L;
  for (i=0;i<(n/32000L);i++)
  k+=fread(ptr+i*32000L,size,(size_t)j,stream);
  j=n%32000L;
  k+=fread(ptr+i*32000L,size,(size_t)j,stream);
  return(k);
}
/******************************************************************************/


/*
This procedure computes and simulates the new eta and ystar
*/
void comp_etaY_grid (int n, int p, int S, int order, int m, int Nknots,
		     double *knots, double *knotsI,  double *G, double *all_beta, double *all_sigma,
		     double *data, double *eta, double *ystar, double ***Pijs, int dens_comp)
{
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

/*printing out Pij
Rprintf("Pij=:");
for (i=0; i<n; i++){
  Rprintf( "i=%d\n", i );
    for (j=0; j<p; j++){
      Rprintf("%2.9f,   ", Pij[i][j]);
    }
   Rprintf( "\n");
}  */

/*printing out Pjs
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


SEXP splb_check(SEXP Rknots, SEXP Rnknots, SEXP RknotsI, SEXP Rorder, SEXP Reta, SEXP Rm,
                SEXP Rallbeta, SEXP Rn, SEXP Rp, SEXP RS, SEXP RG, SEXP Rall_sigma, SEXP Rdata, SEXP Rdens_comp)
{
/* this function coputes Pijs, and selects eta and ystar from the grid and data respectively 
# R check
dyn.load("fcts_Eta_B_win.so")
mu.all.my <- .Call("splb_check", Aknots, as.integer(nknots), AknotsI, as.integer(ord), eta, 
as.integer(m), Rallbeta, as.integer(n), as.integer(p), as.integer(S), grid.eta, Rsigma.vec, data, as.integer(0))
dyn.unload("fcts_Eta_B_win.so")
library(splines)
source("fcts_eta_Bar.R")
eta.R <- comp.etamat.grid(B.ar=B.ar, data, beta, n, p, m, S, grid=grid.eta, sigma.vec=Rsigma.vec, eta.mat=eta)
max(abs(eta.R$prob.arr[,,1]-mu.all.my[[1]]))
max(abs(eta.R$prob.arr[,,2]-mu.all.my[[2]]))
max(abs(eta.R$prob.arr[,,3]-mu.all.my[[3]]))
max(abs(eta.R$prob.arr[,,4]-mu.all.my[[4]]))
max(abs(eta.R$prob.arr[,,5]-mu.all.my[[5]]))
max(abs(eta.R$prob.arr[,,6]-mu.all.my[[6]]))

*/
    int nknots, order, m, n, i, p, nbeta, S;
    double *knots, *knotsI, *eta, *all_beta;
    SEXP ans,  Rystar;
    double *ystar, *G, *all_sigma, *data;
    double ***Pijs;
    int dens_comp;

    nknots=INTEGER(Rnknots)[0];
    S=INTEGER(RS)[0];
    dens_comp=INTEGER(Rdens_comp)[0];

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

    comp_etaY_grid (n, p, S, order, m, nknots, knots, knotsI,  G, all_beta, all_sigma, data, eta, ystar, Pijs, dens_comp);/**/

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


/*
    SEXP Reta_mod;
    PROTECT(ans = allocVector(VECSXP, 2));
    PROTECT(Reta_mod = allocMatrix(REALSXP, n, p));
    PROTECT(Rystar = allocMatrix(REALSXP, S, p));
    for(i = 0; i < n*p; i++) {REAL(Reta_mod)[i] = eta[i];}
    for(i = 0; i < S*p; i++) {REAL(Rystar)[i] = ystar[i];}
    SET_VECTOR_ELT(ans, 0, Reta_mod);
    SET_VECTOR_ELT(ans, 1, Rystar);
    UNPROTECT(3);*/

    
int s1, j;
PROTECT(ans = allocVector(VECSXP, S));
    for (s1=0; s1<S; s1++){
      PROTECT(Rystar = allocMatrix(REALSXP, n, p));
      for(j = 0; j < p; j++) { for (i=0; i<n; i++) {
        REAL(Rystar)[j*n+i]=Pijs[i][j][s1];
        SET_VECTOR_ELT(ans, s1, Rystar);}
      }
      UNPROTECT(1);
    }
    UNPROTECT(1);  


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

/*********************************************************************************/
/*
this procedure computes the overall matrix B (Ball_mat) and the total values of mu_1, ..., mu_p (mu)
*/

void mu_B_fct (int n, int p, int order, int m, int nknots, double *knots, 
	       double *knotsI,  double *all_beta, double *eta, double *Ball_mat, double *mu, FILE *out_allmu)
{

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
    
    for (j=1; j<p; j++){
      if (offsets[j]!=0) {
	     for (l=0; l<order; l++)Ball[m*j+offsets[j]+l-1]=B[j*order+l];
      } else {for (l=1; l<order; l++) Ball[m*j+offsets[j]+l-1]=B[j*order+l];}
    } 
    for (j=0;j<p; j++){
      for (k=0; k<=j; k++){  
      	for (supi=k*m; supi<(k+1)*m; supi++) all_mu[j*(j+1)/2+k]+=Ball[supi] * all_beta[j*(j+1)*m/2+supi];
      }
    }
    for (j=0; j<p; j++){for (k=0;k<j+1; k++){ mu_p[j]=mu_p[j]+all_mu[j*(j+1)/2+k]; }}
    for (k=0; k<p; k++) mu[k*n+i]=mu_p[k];
    for (k=0; k<m*p; k++) Ball_mat[k*n+i]=Ball[k];
   /*if (out_allmu) fsafewrite(all_mu,sizeof(double),(p+1)*p/2,  out_allmu); */
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
/*
this procedure computes the overall matrix B (Ball_mat) and the total values of mu_1, ..., mu_p (mu)
*/

void mu_only_fct (int n, int p, int order, int m, double *all_beta, double *Ball_mat, double *mu)
{

  int i, j,k, supi;
  double *all_mu, *mu_p; 
  
  all_mu=(double *) calloc(p*(p+1)/2, sizeof(double));
  mu_p=(double *) calloc(p, sizeof(double));

  for (i=0; i<n; i++){ 
    memset( all_mu,0, p*(p+1)*sizeof(double)/2);
    memset( mu_p,0, p*sizeof(double));
    for (j=0;j<p; j++){
      for (k=0; k<=j; k++){for (supi=k*m; supi<(k+1)*m; supi++) all_mu[j*(j+1)/2+k]+=Ball_mat[supi] * all_beta[j*(j+1)*m/2+supi];}
    }
    for (j=0; j<p; j++){for (k=0;k<j+1; k++){ mu_p[j]=mu_p[j]+all_mu[j*(j+1)/2+k]; }}
    for (k=0; k<p; k++) mu[k*n+i]=mu_p[k];

  }
  free(all_mu);
  free(mu_p);
}

/******************************************************************************************************/
/* written to check values of mu and B*/

SEXP B_mu_check(SEXP Rknots, SEXP Rnknots, SEXP RknotsI, SEXP Rorder, SEXP Reta, SEXP Rm,
                SEXP Rallbeta, SEXP Rn, SEXP Rp, SEXP Rall_sigma, SEXP Rdata)
{
/* evaluate the non-zero B-spline basis functions (or their derivatives)
 * at xvals.
R check:
mu.R <- cbind(B.ar[,1:m,p]%*%beta[1:m],B.ar[,1:(2*m),p]%*%beta[(m+1):(3*m)],
              B.ar[,1:(3*m),p]%*%beta[(3*m+1):(6*m)])
dyn.load("fcts_Eta_B.so")

mu.all.my <- .Call("B_mu_check", Aknots, as.integer(nknots), AknotsI, as.integer(ord), eta,
as.integer(m), Rallbeta, as.integer(n), as.integer(p), Rsigma.vec, data)

dyn.unload("fcts_Eta_B.so")
mu.R -mu.all.my[[2]]
 */
    int nknots, order, m, n, i, p, nbeta;
    double *knots, *knotsI, *eta, *all_beta;
    SEXP ans, RBall_mat, Rmu;/**/
    double *all_sigma, *data, *Ball_mat, *mu;
    FILE *out_allmu; 
    char  f_allmu_name[17]="allmu_draws_nai";


    nknots=INTEGER(Rnknots)[0];

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

    out_allmu=fopen(f_allmu_name, "wb");
    mu_B_fct(n, p, order,m, nknots,knots, knotsI,  all_beta, eta, Ball_mat, mu, out_allmu);
    fclose(out_allmu);

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

/**************************************************************************/



/* this fct simulate sigma_j*/
double sim_sigma_j(double *yj, double *muj, int n, double asigj, double bsigj)
{
  double sigj, sclsigj, shpsigj;
  int i;
  sclsigj=bsigj; /*1/*/
  for (i=0; i<n; i++) sclsigj=sclsigj+0.5*pow(yj[i]-muj[i],2);
  /*Rprintf("Requiv of bsigj=%f\n", sclsigj);*/
  sclsigj=1/sclsigj;
  shpsigj=asigj+0.5*n;
  sigj=1/rgamma(shpsigj, sclsigj);
  return sigj;
}
/* this fct simulate tau_kj*/
double sim_taukj_naive(double ataukj, double btaukj, double *betakj, int m)
{
  double taukj, scltaukj, shptaukj;
  int i;
  scltaukj=btaukj; /*1/*/
  for (i=0; i<m; i++) scltaukj=scltaukj+0.5*pow(betakj[i],2);
  scltaukj=1/scltaukj;
  shptaukj=ataukj+0.5*m;
  taukj=1/rgamma(shptaukj, scltaukj);
  return taukj; 
}

void getXtX(double *X,int *r,int *c, double *XtX)
/* form X'X (nearly) as efficiently as possible 
   r is number of rows, 
   c is number of columns */
{ double *p0,*p1,*p2,*p3,*p4,x;
  int i,j;
  for (p0=X,i=0;i<*c;i++,p0 += *r) 
  for (p1=X,j=0;j<=i;j++,p1 += *r) {
    for (x=0.0,p2=p0,p3=p1,p4=p0 + *r;p2<p4;p2++,p3++) x += *p2 * *p3;    
    XtX[i + j * *c] = XtX[j + i * *c] = x;
  }
}

void transpose_mat(double *mat, int ncol, int nrow, double *resmat)
{
 int j=0, i=0;
       for(i = 0; i < nrow; i++)    /* column number for the element of resmat */
         for (j = 0; j < ncol; j++)     /* line number for the element of resmat */
           {
              resmat[i* ncol+j]=mat[j* nrow+i];
           }
}

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

void post_betaj_mean_cov(double *Bj_mat, double *tauj, double *yj, double sigj, int j, int m, int n, double *meanj,  double *half_covj) 
/*
* this function computes the peoterior mean and covariance for beta_j
*  meanj is a mj vector and covj,half_covj mj x mj matrices
* 
*/
{
  double *BtB, *Bt, *sigtauj1, *diagcovj, *covj;
  int i,k, mj1, get_inv;
  mj1=m*(j+1);
  BtB=(double *) calloc(mj1*mj1, sizeof(double));
  Bt=(double *) calloc(mj1*n, sizeof(double));
  sigtauj1=(double *) calloc((j+1), sizeof(double));
  diagcovj=(double *) calloc(mj1,sizeof(double));
  getXtX(Bj_mat,&n,&mj1, BtB);
  covj=(double*)calloc(mj1*mj1, sizeof(double));
  /*print_matrix("BtB", mj1,mj1, BtB, mj1);*/
  for (i=0; i<(j+1); i++) sigtauj1[i]=sigj/tauj[i];
  /*print_matrix("sigtauj1", (j+1), 1,sigtauj1, 1);*/
  
  for (i=0; i<(j+1); i++) for (k=0; k<m;k++) BtB[mj1*m*i+mj1*k+m*i+k]=BtB[mj1*m*i+mj1*k+m*i+k]+sigtauj1[i];
  
  get_inv=1;
  /*Rprintf("j+1=%d", j+1);
  print_matrix("BtB", mj1,mj1, BtB, mj1);
  Rprintf("\n");*/
  qr_ldet_inv(BtB,&mj1,covj,&get_inv, diagcovj);/**/
  for (i=0; i<mj1*mj1; i++) half_covj[i]=covj[i];  /*covj is now computed and half_covj is currently equal to covj*/
  free(sigtauj1);
  free(diagcovj);

  
  BtB=(double *) realloc(BtB, mj1*n*sizeof(double));
  transpose_mat(Bj_mat, mj1, n, Bt); /*Bt now contains t(Bj_mat)*/
  matprod(covj, mj1, mj1, Bt,  n, mj1, BtB); /*BtB now contains (BtB)^{-1}%*%t(Bj_mat)*/
  matprod(BtB,  n, mj1, yj, 1, n, meanj); /*meanj now contains (BtB)^{-1}%*%t(Bj_mat)%*%yj */
  /*print_matrix("meanj", n,1, meanj, 1);*/
  free(Bt);
  free(covj);
  
  mroot(half_covj,&mj1,&mj1);/*now half_covj contains (covj)^0.5 but stored by rows*/
  BtB=(double *) realloc(BtB, mj1*mj1*sizeof(double));
  transpose_mat(half_covj, mj1, mj1, BtB);
  for (i=0; i<mj1*mj1; i++) half_covj[i]=BtB[i];/**/
  free(BtB);
}

/****************************************************************************************************************************************************/
/* simulate beta_j given meanj and half_covj*/
void sim_betaj(int j, int m, double *meanj, double *half_covj, double *betaj)
{
  int k, mj1;
  double *hlp;
  mj1=m*(j+1);
  hlp=(double *) calloc(mj1, sizeof(double));
  for (k=0; k<mj1; k++) betaj[k]=rnorm(0, 1);
  matprod(half_covj, mj1, mj1, betaj,  1, mj1, hlp);
  for (k=0; k<mj1; k++) betaj[k]=hlp[k]+meanj[k];
  free(hlp);
}

/****************************************************************************************************************************************************/
/* written to check whether mean and covariance of beta_j posterior are correct */
SEXP betaj_post_check(SEXP Rknots, SEXP Rnknots, SEXP RknotsI, SEXP Rorder, SEXP Reta, SEXP Rm,
                SEXP Rallbeta, SEXP Rn, SEXP Rp, SEXP Rall_sigma, SEXP Rdata, SEXP Rtau)
{
/* this function checks whether we compute correctly parameters for betaj posterior distribution

#R checking
dyn.load("fcts_Eta_B_win.dll")
Rtau <- c(rep(1,p))
mu.all.my <- .Call("betaj_post_check", Aknots, as.integer(nknots), AknotsI, as.integer(ord), eta,
                   as.integer(m), Rallbeta, as.integer(n), as.integer(p), Rsigma.vec, data, Rtau)

dyn.unload("fcts_Eta_B_win.dll")

sigma.vec <- c()
tauj <- c()
for (j in 1:p){
  yj <- as.vector(data[,j])
  Bj <- B.ar[,1:(m*j),j]
  sigmaj<-1
  sigma.vec<-c(sigma.vec,sigmaj)
  tauj <- c(tauj, 1)

  if (length(tauj)>1) T.1=(diag(1/tauj))%x%Im else T.1=Im/tauj
  V.betaj<-solve(sigmaj*T.1+t(Bj)%*%Bj)
  betaj.hat<-V.betaj%*%t(Bj)%*%yj
  cat(max(abs(mu.all.my[[(j-1)*3+1]]-V.betaj)), "\n")
  cat( max(abs(mu.all.my[[(j-1)*3+3]]-betaj.hat)), "\n")
 }
*/
    int nknots, order, m, n, i, j, p, nbeta, mj1;
    double *knots, *knotsI, *eta, *all_beta;
    SEXP ans, Rcov, Rhalfcov, Rmean;
    double *all_sigma, *data, *Ball_mat, *Bj_mat, *tauj;
    double sigj,*mu, *meanj,*half_covj, *yj; 
    FILE *out_allmu; 
    char  f_allmu_name[17]="allmu_draws_nai";



    nknots=INTEGER(Rnknots)[0];
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
    
    out_allmu=fopen(f_allmu_name, "wb");    
    mu_B_fct(n, p, order, m, nknots,knots, knotsI,  all_beta, eta, Ball_mat, mu, out_allmu);
    fclose(out_allmu);
    
    Bj_mat=(double*)calloc(m*n, sizeof(double));
    tauj=(double*)calloc(1, sizeof(double));
    yj=(double*)calloc(n, sizeof(double));
    meanj=(double*)calloc(m, sizeof(double));
    half_covj=(double*)calloc(m*m, sizeof(double));
    PROTECT(ans = allocVector(VECSXP, 9));
    
    for (j=0; j<p; j++) {
      mj1=m*(j+1);
      Bj_mat=(double*)realloc(Bj_mat, mj1*n*sizeof(double));
      for (i=0; i<mj1*n; i++) Bj_mat[i]=Ball_mat[i];
      tauj=(double*)realloc(tauj, (j+1)*sizeof(double));
      for (i=0; i<(j+1); i++) tauj[i]=REAL(Rtau)[i];
      for (i=0; i<n; i++) yj[i]=data[n*j+i];
      meanj=(double*)realloc(meanj, mj1*sizeof(double));
      half_covj=(double*)realloc(half_covj, mj1*mj1*sizeof(double));
      sigj=all_sigma[j];
      post_betaj_mean_cov(Bj_mat, tauj, yj,  sigj,  j,  m,  n, meanj, half_covj);
      PROTECT(Rcov = allocMatrix(REALSXP, mj1, mj1));
      PROTECT(Rhalfcov = allocMatrix(REALSXP, mj1, mj1));
      PROTECT(Rmean = allocVector(REALSXP, mj1));
      for(i = 0; i < mj1*mj1; i++) {REAL(Rhalfcov)[i] = half_covj[i];}
      for(i = 0; i < mj1; i++) REAL(Rmean)[i] = meanj[i];
      SET_VECTOR_ELT(ans, j*3, Rcov);
      SET_VECTOR_ELT(ans, j*3+1, Rhalfcov);
      SET_VECTOR_ELT(ans, j*3+2, Rmean);      

      UNPROTECT(3);
    }
UNPROTECT(1);   

    free(Bj_mat);
    free(tauj);
    free(yj);
    free(meanj);
    free(half_covj);
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
/* Gibbs sampler for BALVM with naive prior */

void one_step_gibbs_naive(int m, int n, int p, double *Ball_mat, double *data, double *mu, double a_sigma, double b_sigma, double a_tau, double b_tau,double *all_beta, double *sigma_vec )
{
 int j,i,s, mj1, k;
 double sig_j, tau_kj;
 double *Bj_mat, *betaj,*betakj, *tauj, *mu_betaj, *half_covj, *yj, *muj;
 
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
    sig_j=sim_sigma_j(yj, muj, n, a_sigma, b_sigma);
    sigma_vec[j]=sig_j;
    Rprintf("j=%d, sig_j=%f\n", j, sig_j);


    for (k=0; k<j+1; k++)
    {
      for (s=0; s<m; s++) betakj[s]=betaj[k*m+s];
      tau_kj=sim_taukj_naive(a_tau, b_tau, betakj, m);
      tauj[k]=tau_kj;
    }
  
  post_betaj_mean_cov(Bj_mat, tauj, yj, sig_j, j, m, n, mu_betaj, half_covj);
  sim_betaj(j, m, mu_betaj, half_covj, betaj);
  for (s=0; s<mj1; s++) all_beta[j*(j+1)*m/2+s]=betaj[s];
  free(Bj_mat);
  free(betaj);
  free(tauj);
  free(mu_betaj);
  free(half_covj);
  }
 free(betakj);
 free(yj);
 free(muj);

}


/*******************************************************************************************************/
/* Gibbs sampler for BALVM with naive prior */

SEXP gibbs_BALVM_naive(SEXP Rtot_draws, SEXP Rknots, SEXP RNknots, SEXP RknotsI, SEXP Rorder,
          SEXP Rn, SEXP Rm, SEXP Rp, SEXP RGrid_length, SEXP Rgrid, SEXP Rdata, SEXP Rab_sigma, SEXP Rab_tau, SEXP Reta_start,SEXP Rallbeta_start, SEXP Rdens_comp, SEXP Rdrws_step)
{
 int i, j,  m, p, n, order, grid_length;
 int pn, t, tot_draws, Nknots;
 int dens_comp, drws_step, it_step, eff_draws;
 double *data, *ystar, *knots, *knotsI;
 double *eta, *all_beta;
 double *mu, *Ball_mat;
 double a_sigma, b_sigma, a_tau, b_tau;
 double *sigma_vec;
 double *grid;
 char f_mu_name[13]="mu_draws_nai", f_sig_name[16]="sigma_draws_nai", 
      f_mu_spec_name[18]="mu_spec_draws_nai", f_sig_spec_name[21]="sigma_spec_draws_nai",
      f_eta_name[13]="eta_draws_nai", f_beta_name[15]="beta_draws_nai", f_allmu_name[17]="allmu_draws_nai";
 
 double *eta_spec, *allbeta_spec, *Ball_mat_spec, *mu_spec;
 
 FILE *out_mu, *out_sigma, *out_mu_spec, 
      *out_sigma_spec,*in_mu, *in_sigma, 
      *in_mu_spec, *in_sigma_spec, 
      *out_eta,  *in_eta, *out_beta,  *in_beta, *out_allmu , *out_allmu_spec;
 
 SEXP ans, Rmu, Rsigma, Rsigma_spec, Reta, Rbeta;
 
 Nknots=INTEGER(RNknots)[0];
 drws_step=INTEGER(Rdrws_step)[0];
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
 dens_comp=INTEGER(Rdens_comp)[0];

 double ***Pijs;
 Pijs=DM3D(n,p,grid_length);

 mu=(double *) malloc(pn*sizeof(double));
 Ball_mat=(double *) calloc(n*m*p,sizeof(double));
 sigma_vec=(double *) malloc(p*sizeof(double));
 data=(double *) malloc(pn*sizeof(double));
 eta=(double *) malloc(pn*sizeof(double));
 ystar=(double *) malloc(grid_length*p*sizeof(double));
 
 for (i=0; i<pn; i++) {data[i]=REAL(Rdata)[i]; eta[i]=REAL(Reta_start)[i];}
 
 grid=(double *) malloc(grid_length*sizeof(double));
 for (i=0; i<grid_length; i++) grid[i]=REAL(Rgrid)[i];
 
 all_beta=(double *) malloc(p*(p+1)*m*sizeof(double)/2);
 for (i=0; i<p*(p+1)*m/2; i++) all_beta[i]=REAL(Rallbeta_start)[i];

 knots=(double *) malloc(Nknots*sizeof(double));
 for (i=0; i<Nknots; i++) knots[i]=REAL(Rknots)[i];
 knotsI=(double *) malloc((Nknots-1)*sizeof(double));
 for (i=0; i<(Nknots-1); i++) knotsI[i]=REAL(RknotsI)[i];
 
 /*out_allmu=fopen(f_allmu_name, "wb");*/
 mu_B_fct (n, p, order, m, Nknots, knots, knotsI,  all_beta, eta, Ball_mat, mu, out_allmu);
  
if (dens_comp==1)
{
 mu_spec=(double *) malloc(p*grid_length*sizeof(double));
 Ball_mat_spec=(double *) calloc(grid_length*m*p,sizeof(double));
 eta_spec=(double *) malloc(p*grid_length*sizeof(double));
 for (j=0; j<p; j++) {for (i=0; i<grid_length; i++) eta_spec[j*grid_length+i]=grid[i];}
 allbeta_spec=(double *) malloc(p*(p+1)*m*sizeof(double)/2);
 for (i=0; i<p*(p+1)*m/2; i++) allbeta_spec[i]=all_beta[i];
 mu_B_fct (grid_length, p, order, m, Nknots, knots, knotsI,  allbeta_spec, eta_spec, Ball_mat_spec, mu_spec,out_allmu_spec);
 out_mu_spec=fopen(f_mu_spec_name, "wb");
 out_sigma_spec=fopen(f_sig_spec_name, "wb");
 /*print_matrix("mu_spec", grid_length, p, mu_spec, grid_length);*/
}

 
 out_mu=fopen(f_mu_name, "wb");
 out_sigma=fopen(f_sig_name, "wb");
 out_eta=fopen(f_eta_name, "wb");
 out_beta=fopen(f_beta_name, "wb");
 
 GetRNGstate();
 /*starting draws*/
 it_step=0;

for (eff_draws=0; eff_draws<tot_draws; )
 {
 
 one_step_gibbs_naive(m, n, p, Ball_mat, data, mu, a_sigma, b_sigma, a_tau, b_tau, all_beta, sigma_vec );
 /*print_matrix("sigma_vec, non density", p,1,sigma_vec,1);*/
 comp_etaY_grid(n, p, grid_length, order, m, Nknots, knots, knotsI, grid, all_beta, sigma_vec, data, eta, ystar, Pijs, dens_comp);
 mu_B_fct (n, p, order, m, Nknots, knots, knotsI,  all_beta, eta, Ball_mat, mu, out_allmu); 
 it_step=it_step+1;
 if (it_step==drws_step)
 {
  fsafewrite(mu, sizeof(double), pn, out_mu);
  print_matrix("sigma avant fsafewrite",p, 1,sigma_vec, p);/**/
  fsafewrite(sigma_vec, sizeof(double),p, out_sigma);
  print_matrix("sigma apres fsafewrite",p, 1,sigma_vec, p);/**/
  fsafewrite(eta, sizeof(double),pn, out_eta);
  fsafewrite(all_beta, sizeof(double),p*(p+1)*m/2, out_beta);
   if (dens_comp==0)
   {
     eff_draws=eff_draws+1;
     it_step=0;     
   }
 }
 
 if (dens_comp==1)
  {
  one_step_gibbs_naive(m, grid_length, p, Ball_mat_spec, ystar, mu_spec,  a_sigma, b_sigma, a_tau, b_tau, allbeta_spec, sigma_vec );
  print_matrix("sigma_vec, density", p,1,sigma_vec,1);
  mu_B_fct (grid_length, p, order, m, Nknots, knots, knotsI,  allbeta_spec, eta_spec, Ball_mat_spec, mu_spec,out_allmu_spec); 
  /*mu_only_fct (grid_length, p, order, m, allbeta_spec, Ball_mat_spec, mu_spec);
  print_matrix("mu_spec", grid_length, p, mu_spec, grid_length);*/
 
  if (it_step==drws_step)
  {
    fsafewrite(mu_spec, sizeof(double), p*grid_length, out_mu_spec);
    fsafewrite(sigma_vec, sizeof(double),p, out_sigma_spec);
    eff_draws=eff_draws+1;
    it_step=0; 
  }
 }
 Rprintf("it_step=%d,  eff_draws=%d\n",it_step, eff_draws);
 } 
 PutRNGstate();

 /*fclose(out_allmu);*/
 fclose(out_mu);
 fclose(out_sigma);
 fclose(out_eta);
 fclose(out_beta);
 
 if(dens_comp==1)
 {
  fclose(out_mu_spec);
  fclose(out_sigma_spec);
 }

 free(data);
 free(grid);
 free(knots);
 free(knotsI);
 free(Ball_mat);
 free(ystar);
 Dfree3d(Pijs);
 


 if (dens_comp==1)
  {
  free(eta_spec);
  free(allbeta_spec);
  free(Ball_mat_spec);
  }
 
 in_mu=fopen(f_mu_name, "rb");
 in_sigma=fopen(f_sig_name, "rb");
 in_eta=fopen(f_eta_name, "rb");
 in_beta=fopen(f_beta_name, "rb");
 
 if (dens_comp==1)
  {
    in_mu_spec=fopen(f_mu_spec_name, "rb");
    in_sigma_spec=fopen(f_sig_spec_name, "rb");
  }

  if (dens_comp==1)
  {
   PROTECT(ans = allocVector(VECSXP, 3*tot_draws+2));
   PROTECT(Rsigma = allocMatrix(REALSXP,  p, tot_draws));
   PROTECT(Rsigma_spec = allocMatrix(REALSXP,  p, tot_draws));
  } else
    
  {
   PROTECT(ans = allocVector(VECSXP, 3*tot_draws+1));
   PROTECT(Rsigma = allocMatrix(REALSXP,  p, tot_draws));
  }
 for (t=0; t<tot_draws; t++) {
    fsaferead(mu, sizeof(double), pn, in_mu);
    fsaferead(sigma_vec, sizeof(double), p, in_sigma);
    fsaferead(eta, sizeof(double), pn, in_eta);
    fsaferead(all_beta, sizeof(double), p*(p+1)*m/2, in_beta);
    PROTECT(Rmu = allocMatrix(REALSXP, n, p));
    PROTECT(Reta = allocMatrix(REALSXP, n, p));
    PROTECT(Rbeta = allocVector(REALSXP, p*(p+1)*m/2));
    for(i = 0; i < pn; i++) {REAL(Rmu)[i] = mu[i]; REAL(Reta)[i] = eta[i];}
    for(i = 0; i < p*(p+1)*m/2; i++) {REAL(Rbeta)[i] = all_beta[i];}
    for(i = 0; i < p; i++) REAL(Rsigma)[t*p+i] = sigma_vec[i];
    SET_VECTOR_ELT(ans, t, Rmu);
    SET_VECTOR_ELT(ans, tot_draws+t, Reta);
    SET_VECTOR_ELT(ans, 2*tot_draws+t, Rbeta);
    UNPROTECT(3);
    
    if (dens_comp==1)
    {
      fsaferead(mu_spec, sizeof(double), p*grid_length, in_mu_spec);
      fsaferead(sigma_vec, sizeof(double), p, in_sigma_spec);
      PROTECT(Rmu = allocMatrix(REALSXP, grid_length, p));
      for(i = 0; i < p*grid_length; i++) REAL(Rmu)[i] = mu_spec[i];
      for(i = 0; i < p; i++) REAL(Rsigma_spec)[t*p+i] = sigma_vec[i];
      SET_VECTOR_ELT(ans, 2*tot_draws+t, Rmu);
      UNPROTECT(1);
    }
 }
 if (dens_comp==1) {
  SET_VECTOR_ELT(ans, 3*tot_draws, Rsigma);
  SET_VECTOR_ELT(ans, 3*tot_draws+1, Rsigma_spec);
  UNPROTECT(3);
  fclose(in_mu);
  fclose(in_sigma);
  fclose(in_mu_spec);
  fclose(in_sigma_spec);
  free(mu);
  free(mu_spec);
  free(sigma_vec);
 } else {
   SET_VECTOR_ELT(ans, 3*tot_draws, Rsigma);
   UNPROTECT(2);
    fclose(in_mu);
    fclose(in_sigma);
    fclose(in_eta);
    fclose(in_beta);
    free(mu);
    free(sigma_vec);
    free(eta);
    free(all_beta);
 }
 return (ans);
}

/***************************************************************************************************/
void Kkj_mat(int m, double *Kkj)
{
  int i;
 Kkj[0]=1.0;
 Kkj[1]=-1.0;
 Kkj[m*m-2]=-1.0;
 Kkj[m*m-1]=1.0;
 for (i=1; i<(m-1); i++)
 {
   Kkj[i*m+i-1]=-1.0;
   Kkj[i*m+i]=2.0;
   Kkj[i*m+i+1]=-1.0;
 }
}

/***************************************************************************************************/
void scalar_mult_Kkj_mat(double sc, int m, double *Kkj)
{
  int i;
 Kkj[0]=Kkj[0]*sc;
 Kkj[1]=Kkj[1]*sc;
 Kkj[m*m-2]=Kkj[m*m-2]*sc;
 Kkj[m*m-1]=Kkj[m*m-1]*sc;
 for (i=1; i<(m-1); i++)
 {
   Kkj[i*m+i-1]=Kkj[i*m+i-1]*sc;
   Kkj[i*m+i]=Kkj[i*m+i]*sc;
   Kkj[i*m+i+1]=Kkj[i*m+i+1]*sc;
 }
}

/**************************************************************************************************************************************************************/

void Psplines_betakj_mean_cov(double *Bkj_mat, double taukj, double *muj_k, double *yj, double sigj, int m, int n, double *meanj, double *half_covj)
/*
* this function computes the peoterior mean and covariance for beta_j
*  meanj is a mj vector and covj,half_covj mj x mj matrices
* 
*/
{
  double *BtB, *Bt, *Kkj, *diagcovj, *covj;
  int i,get_inv;
  BtB=(double *) calloc(m*m, sizeof(double));
  Bt=(double *) calloc(m*n, sizeof(double));
  Kkj=(double *) calloc(m*m, sizeof(double));
  diagcovj=(double *) calloc(m,sizeof(double));
  covj=(double*)calloc(m*m, sizeof(double));
  Kkj_mat(m, Kkj);
  scalar_mult_Kkj_mat(1/taukj, m, Kkj);
  getXtX(Bkj_mat,&n,&m, BtB);
  for (i=0; i<m*m; i++) BtB[i]=BtB[i]+Kkj[i];
  
  get_inv=1;
  /*Rprintf("j+1=%d", j+1);
  print_matrix("BtB", mj1,mj1, BtB, mj1);
  Rprintf("\n");*/
  qr_ldet_inv(BtB,&m,covj,&get_inv, diagcovj);/**/
  for (i=0; i<m*m; i++) half_covj[i]=covj[i];  /*covj is now computed and half_covj is currently equal to covj*/
  
  BtB=(double *) realloc(BtB, m*n*sizeof(double));
  transpose_mat(Bkj_mat, m, n, Bt); /*Bt now contains t(Bj_mat)*/
  matprod(covj, m, m, Bt,  n, m, BtB); /*BtB now contains (BtB/sigj+Kkj/taukj)^{-1}%*%t(Bj_mat)*/
  diagcovj=(double *) realloc(diagcovj, n*sizeof(double));
  for (i=0; i<n; i++) diagcovj[i]= (yj[i]-muj_k[i])/sigj;
  /*print_matrix("yj-muj_k=", n,1, diagcovj, 1);*/
  matprod(BtB,  n, m, diagcovj, 1, n, meanj); /*meanj now contains (BtB)^{-1}%*%t(Bj_mat)%*%(yj-muj_k) */
  /*print_matrix("meanj", n,1, meanj, 1);*/
  free(Bt);
  
  mroot(half_covj,&m,&m);/*now half_covj contains (covj)^0.5 but stored by rows*/
  BtB=(double *) realloc(BtB, m*m*sizeof(double));
  transpose_mat(half_covj, m, m, BtB);
  for (i=0; i<m*m; i++) half_covj[i]=BtB[i];/**/
  free(BtB);
  free(diagcovj);
  free(covj);
  free(Kkj);
}

double trBtAB(double *A,double *B,int *n,int*m) 
/* form tr(B'AB) where A is n by n and B is n by m, m < n,
   basic point is that this is sum_ijk A_ik B_ij B_kj
 */
{ double tr=0.0,x,*p,*p1,*p2;
  int j,k;
  for (j=0;j<*m;j++)
    for (k=0;k<*n;k++) {
      p = A + *n * k;p2 = p + *n;
      p1 = B + *n * j;
      x = B[k + j * *n];
      for (;p<p2;p++,p1++) tr += *p * *p1 * x;   
    }
  return(tr);
}
/*******************************************************************************************************/
/* this fct simulate tau_kj*/

double sim_taukj_psplines(double ataukj, double btaukj, double *betakj, double *Kkj, int m)
{
  double taukj, scltaukj, shptaukj;
  int one=1;
  scltaukj=btaukj; /*1/*/
  scltaukj=scltaukj+0.5*trBtAB(Kkj, betakj, &m, &one);
  /*Rprintf("Requiv of btaukj=%f\n", scltaukj);*/
  scltaukj=1/scltaukj;
  shptaukj=ataukj+0.5*(m-1);
  taukj=1/rgamma(shptaukj, scltaukj);
  return taukj; 
}

 /*******************************************************************************************************/
 
 void sim_betakj_psplines(int m, double *mu_betakj, double *half_covkj, double *betakj)
{
  int k;
  double *hlp;
  hlp=(double *) calloc(m, sizeof(double));
  
  for (k=0; k<m; k++) betakj[k]=rnorm(0, 1);
  
  matprod(half_covkj, m, m, betakj, 1, m, hlp);
  
  for (k=0; k<m; k++) betakj[k]=hlp[k]+mu_betakj[k];
  
  free(hlp);
}
/*******************************************************************************************************/
/* one step of Gibbs sampler for BALVM with P-splines prior */

void psplines_one_step(int m, int n, int p, double *Ball_mat, double *data, double *mu, double a_sigma, double b_sigma, double a_tau, double b_tau,double *all_beta, double *sigma_vec )
{
int i,j,k,s;
double *yj, *muj, *Bkj_mat, *half_covkj, *betakj, *mu_betakj, *muj_k, *Kkj;
double sig_j, tau_kj;

 Kkj=(double*)calloc(m*m, sizeof(double));
 Kkj_mat(m, Kkj);
 mu_betakj=(double*)calloc(m, sizeof(double));
 betakj=(double*)malloc(m*sizeof(double));
 half_covkj=(double*)calloc(m*m, sizeof(double));
 yj=(double*)malloc(n*sizeof(double));
 muj=(double*)malloc(n*sizeof(double));
 muj_k=(double *) malloc(n*sizeof(double));
 Bkj_mat=(double *) calloc( n*m, sizeof(double));

 for (j=0; j<p; j++)
  {
    for (i=0; i<n; i++) yj[i]=data[n*j+i];
    for (i=0; i<n; i++) muj[i]=mu[j*n+i];
    sig_j=sim_sigma_j(yj, muj, n, a_sigma, b_sigma);
    sigma_vec[j]=sig_j;
    for (k=0; k<j+1; k++)
    {
      for (s=0; s<m; s++) betakj[s]=all_beta[j*(j+1)*m/2+k*m+s];
      tau_kj= sim_taukj_psplines(a_tau, b_tau, betakj, Kkj, m);
      for (s=0; s<m*n; s++) Bkj_mat[s]=Ball_mat[m*n*k+s];
      matprod(Bkj_mat, m, n, betakj,  1, m, muj_k);
      for (s=0; s<n; s++) muj_k[s]=muj[s]-muj_k[s];
      memset( half_covkj, 0, m*m*sizeof(double));
      memset(mu_betakj, 0, m*sizeof(double));
      Psplines_betakj_mean_cov(Bkj_mat, tau_kj, muj_k, yj, sig_j, m, n, mu_betakj, half_covkj);
      memset(betakj, 0, m*sizeof(double));
      sim_betakj_psplines(m, mu_betakj, half_covkj, betakj);
      for (s=0; s<m; s++) all_beta[j*(j+1)*m/2+k*m+s]=betakj[s];

    }

  }
    free(yj); free(muj);
    free(Bkj_mat);  free(muj_k);
    free(mu_betakj);
    free(betakj);
    free(half_covkj);
   }
/*******************************************************************************************************/
/* Gibbs sampler for BALVM with P-splines prior */

SEXP gibbs_BALVM_psplines(SEXP Rtot_draws, SEXP Rknots, SEXP RNknots, SEXP RknotsI, SEXP Rorder,
          SEXP Rn, SEXP Rm, SEXP Rp, SEXP RGrid_length, SEXP Rgrid, SEXP Rdata, SEXP Rab_sigma, SEXP Rab_tau, SEXP Reta_start, SEXP Rallbeta_start, SEXP Rdens_comp, SEXP Rdrws_step)
{
 int i, j, m, p, n, order, grid_length;
 int pn, t, tot_draws, Nknots;
 int dens_comp, drws_step, it_step, eff_draws;
 double *data, *ystar, *knots, *knotsI;
 double *eta, *all_beta;
 double *mu, *Ball_mat;
 double a_sigma, b_sigma, a_tau, b_tau;
 double *sigma_vec;
 double *grid;
 char f_mu_name[13]="mu_draws_psp", f_sig_name[16]="sigma_draws_psp", 
      f_mu_spec_name[18]="mu_spec_draws_psp", f_sig_spec_name[21]="sigma_spec_draws_psp",
      f_beta_name[15]="beta_draws_psp", f_eta_name[14]="eta_draws_psp",f_allmu_name[17]="allmu_draws_psp";
      
 double *eta_spec, *allbeta_spec, *Ball_mat_spec, *mu_spec;
 FILE *out_mu, *out_sigma, 
      *out_mu_spec,*out_eta, 
      *out_sigma_spec,*out_beta, 
      *in_mu, *in_sigma, *in_mu_spec, 
      *in_sigma_spec, *in_eta, *in_beta, *out_allmu, *out_allmu_spec;
      
 SEXP ans, Rmu, Rsigma, Rsigma_spec, Reta, Rbeta;
 
 Nknots=INTEGER(RNknots)[0];
 drws_step=INTEGER(Rdrws_step)[0];
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
 dens_comp=INTEGER(Rdens_comp)[0];
 
 data=(double *) malloc(pn*sizeof(double));
 eta=(double *) malloc(pn*sizeof(double));
 
 for (i=0; i<pn; i++) {data[i]=REAL(Rdata)[i]; eta[i]=REAL(Reta_start)[i];}
 
 grid=(double *) malloc(grid_length*sizeof(double));
 for (i=0; i<grid_length; i++) grid[i]=REAL(Rgrid)[i];
 
 all_beta=(double *) malloc(p*(p+1)*m*sizeof(double)/2);
 for (i=0; i<p*(p+1)*m/2; i++) all_beta[i]=REAL(Rallbeta_start)[i];

 knots=(double *) malloc(Nknots*sizeof(double));
 for (i=0; i<Nknots; i++) knots[i]=REAL(Rknots)[i];
 knotsI=(double *) malloc((Nknots-1)*sizeof(double));
 for (i=0; i<(Nknots-1); i++) knotsI[i]=REAL(RknotsI)[i];
 
 mu=(double *) malloc(pn*sizeof(double));
 Ball_mat=(double *) calloc(n*m*p,sizeof(double));
 /*print_matrix("all_beta=",p*(p+1)*m/2, 1,all_beta, 1 );
 print_matrix("knots=",Nknots, 1,knots, 1 );
 print_matrix("knotsI=",Nknots-1, 1,knotsI, 1 );*/
 out_allmu=fopen(f_allmu_name, "wb");
 mu_B_fct (n, p, order, m, Nknots, knots, knotsI,  all_beta, eta, Ball_mat, mu, out_allmu);
 /*print_matrix("mu=",n*p,1, mu,1 );*/
 
 sigma_vec=(double *) malloc(p*sizeof(double));
 ystar=(double *) malloc(grid_length*p*sizeof(double));
 
 double ***Pijs;
 Pijs=DM3D(n,p,grid_length);

  if (dens_comp==1)
  {
   mu_spec=(double *) calloc(p*grid_length, sizeof(double));
   Ball_mat_spec=(double *) calloc(grid_length*m*p,sizeof(double));
   eta_spec=(double *) malloc(p*grid_length*sizeof(double));
   for (j=0; j<p; j++) {for (i=0; i<grid_length; i++) eta_spec[j*grid_length+i]=grid[i];}
   allbeta_spec=(double *) malloc(p*(p+1)*m*sizeof(double)/2);
   for (i=0; i<p*(p+1)*m/2; i++) allbeta_spec[i]=all_beta[i];
   mu_B_fct (grid_length, p, order, m, Nknots, knots, knotsI,  allbeta_spec, eta_spec, Ball_mat_spec, mu_spec, out_allmu_spec);
   out_mu_spec=fopen(f_mu_spec_name, "wb");
   out_sigma_spec=fopen(f_sig_spec_name, "wb");
  }

 
 out_mu=fopen(f_mu_name, "wb");
 out_sigma=fopen(f_sig_name, "wb");
 out_eta=fopen(f_eta_name, "wb");
 out_beta=fopen(f_beta_name, "wb");
 GetRNGstate();


 it_step=0;
 for (eff_draws=0; eff_draws<tot_draws; )
 {
   psplines_one_step(m, n, p, Ball_mat,data, mu, a_sigma, b_sigma, a_tau, b_tau,all_beta, sigma_vec );
   comp_etaY_grid(n, p, grid_length, order, m, Nknots, knots, knotsI, grid, all_beta, sigma_vec, data, eta, ystar, Pijs, dens_comp);
   mu_B_fct (n, p, order, m, Nknots, knots, knotsI,  all_beta, eta, Ball_mat, mu, out_allmu);
   it_step=it_step+1;
   if (it_step==drws_step)
   {
    fsafewrite(mu, sizeof(double), pn, out_mu);
    fsafewrite(sigma_vec, sizeof(double),p, out_sigma);
    fsafewrite(eta, sizeof(double),pn, out_eta);
    fsafewrite(all_beta, sizeof(double),p*(p+1)*m/2, out_beta);
    if (dens_comp==0)
     {
       it_step=0;
       eff_draws=eff_draws+1;
     }
   }

  if (dens_comp==1)
  {
   psplines_one_step(m, grid_length, p, Ball_mat_spec, ystar, mu_spec, a_sigma, b_sigma, a_tau, b_tau,allbeta_spec, sigma_vec );
   mu_B_fct (grid_length, p, order, m, Nknots, knots, knotsI,  allbeta_spec, eta_spec, Ball_mat_spec, mu_spec, out_allmu_spec);

  if (it_step==drws_step)
 {
  fsafewrite(mu_spec, sizeof(double), p*grid_length, out_mu_spec);
  fsafewrite(sigma_vec, sizeof(double),p, out_sigma_spec);
  it_step=0;
  eff_draws=eff_draws+1;
 }
 }
Rprintf("it_step=%d,  eff_draws=%d\n",it_step, eff_draws);
 } 

 PutRNGstate();
 fclose(out_allmu);
 fclose(out_mu);
 fclose(out_sigma);
 fclose(out_eta);
 fclose(out_beta);

 if(dens_comp==1)
 {
  fclose(out_mu_spec);
  fclose(out_sigma_spec);
 }

 free(data);
 free(grid);
 free(knots);
 free(knotsI);
 free(Ball_mat);
 Dfree3d(Pijs);
 free(ystar);


 if (dens_comp==1)
  {
  free(eta_spec);
  free(allbeta_spec);
  free(Ball_mat_spec);
  }
 
 in_mu=fopen(f_mu_name, "rb");
 in_sigma=fopen(f_sig_name, "rb");
 in_eta=fopen(f_eta_name, "rb");
 in_beta=fopen(f_beta_name, "rb");
 
 if (dens_comp==1)
  {
    in_mu_spec=fopen(f_mu_spec_name, "rb");
    in_sigma_spec=fopen(f_sig_spec_name, "rb");

  }

  if (dens_comp==1)
  {
   PROTECT(ans = allocVector(VECSXP, 3*tot_draws+2));
   PROTECT(Rsigma = allocMatrix(REALSXP,  p, tot_draws));
   PROTECT(Rsigma_spec = allocMatrix(REALSXP,  p, tot_draws));
  } else
  {
   PROTECT(ans = allocVector(VECSXP, 3*tot_draws+1));
   PROTECT(Rsigma = allocMatrix(REALSXP,  p, tot_draws));
  }
 for (t=0; t<tot_draws; t++) {
    fsaferead(mu, sizeof(double), pn, in_mu);
    fsaferead(sigma_vec, sizeof(double), p, in_sigma);
    fsaferead(eta, sizeof(double), pn, in_eta);
    fsaferead(all_beta, sizeof(double), p*(p+1)*m/2, in_beta);
    PROTECT(Rmu = allocMatrix(REALSXP, n, p));
    PROTECT(Reta = allocMatrix(REALSXP, n, p));
    PROTECT(Rbeta = allocVector(REALSXP, p*(p+1)*m/2));
    for(i = 0; i < pn; i++) {REAL(Rmu)[i] = mu[i]; REAL(Reta)[i] = eta[i];}
    for(i = 0; i < p; i++) REAL(Rsigma)[t*p+i] = sigma_vec[i];
    for(i = 0; i < p*(p+1)*m/2; i++) {REAL(Rbeta)[i] = all_beta[i];}
    SET_VECTOR_ELT(ans, t, Rmu);
    SET_VECTOR_ELT(ans, tot_draws+t, Reta);
    SET_VECTOR_ELT(ans, 2*tot_draws+t, Rbeta);
    UNPROTECT(3);
    if (dens_comp==1)
    {
      fsaferead(mu_spec, sizeof(double), p*grid_length, in_mu_spec);
      fsaferead(sigma_vec, sizeof(double), p, in_sigma_spec);
      PROTECT(Rmu = allocMatrix(REALSXP, grid_length, p));
      for(i = 0; i < p*grid_length; i++) REAL(Rmu)[i] = mu_spec[i];
      for(i = 0; i < p; i++) REAL(Rsigma_spec)[t*p+i] = sigma_vec[i];
      SET_VECTOR_ELT(ans, 2*tot_draws+t, Rmu);
      UNPROTECT(1);
    }
 }
 if (dens_comp==1) {
  SET_VECTOR_ELT(ans, 3*tot_draws, Rsigma);
  SET_VECTOR_ELT(ans, 3*tot_draws+1, Rsigma_spec);
  UNPROTECT(3);
  fclose(in_mu);
  fclose(in_sigma);
  fclose(in_mu_spec);
  fclose(in_sigma_spec);
  free(mu);
  free(mu_spec);
  free(sigma_vec);
 } else {
   SET_VECTOR_ELT(ans, 3*tot_draws, Rsigma);
   UNPROTECT(2);
    fclose(in_mu);
    fclose(in_sigma);
    fclose(in_eta);
    fclose(in_beta);
    free(mu);
    free(sigma_vec);
    free(eta);
    free(all_beta);
 }
 return (ans);
}

/**************************************************************************/

/* this fct simulate sigma_j in GDP algorithm*/
double sigmaj_gdp_update(double *yj, double *muj, int n, int j, int m, double *tauj, double *betaj)
{
  double sclsigj, shpsigj,sigmaj, rt;
  int i,k;
  shpsigj=(n+(j+1)*m)/2; /*1/*/
  sclsigj=0;
  if (j==2){
  print_matrix("yj from sigmaj_gdp_update", n ,1, yj,n);
  print_matrix("muj from sigmaj_gdp_update", n ,1, muj,n);}
  /*print_matrix("tauj from sigmaj_gdp_update", (j+1)*(m-1) ,1, tauj,(j+1)*(m-1));
  print_matrix("betaj from sigmaj_gdp_update", (j+1)*m ,1, betaj,(j+1)*m);}*/
  
  for (i=0; i<n; i++) sclsigj=sclsigj+0.5*pow(yj[i]-muj[i],2);
  Rprintf("sclsigj due to (yj-muj)^2=%f\n", sclsigj);
  
  for (k=0; k<(j+1); k++)
  {
    for (i=0; i<m-1; i++) {sclsigj=sclsigj+0.5*pow(betaj[k*m+i+1]-betaj[k*m+i],2)/tauj[k*(m-1)+i];}
  }
  rt=1/sclsigj;
  sigmaj=1/rgamma(shpsigj, rt);
  if (sigmaj>HIGH) sigmaj=HIGH;
  return sigmaj;
}

/**************************************************************************/

double alphaj_gdp_update(double sigmaj, double nuj,double *betaj, int j, int m, int grid_length)
{
  double *a;
  int *perm, *num_grid;
  double prod_beta, suma, alpha;
  double *grid0;
  int jm, i, s;
  grid0=(double *) calloc(grid_length, sizeof(double));
  
  for (i=0;i<grid_length; i++) grid0[i]=0.00001+i*(1-0.00002)/grid_length;

  jm=j*m;
  prod_beta=1.0;
  suma=0;
  for (i=0; i<j; i++) {
    for (s=1; s<m; s++) prod_beta=prod_beta*(1+fabs(betaj[i*m+s]-betaj[i*m+s-1])/(sqrt(sigmaj)*nuj));
  }
  a=(double *) calloc(grid_length, sizeof(double));
  for (i=0; i<grid_length; i++)
  {
   a[i]=pow((1-grid0[i])/grid0[i], jm-j)*pow(prod_beta, -1/grid0[i]);
   /*Rprintf("ISNAN(a[i])=%d\n", ISNAN(a[i]));*/
   if(!R_FINITE(a[i])||ISNAN(a[i])||ISNA(a[i])) a[i]=0;
   suma=suma+a[i];
  }
  if (suma>0) {for (i=0; i<grid_length; i++) a[i]=a[i]/suma;} else {for (i=0; i<grid_length; i++) a[i]=1/grid_length;}
  /*print_matrix("a mat from fct alphaj_gdp_update", grid_length,1,a,1 );*/
  perm=(int*) calloc(grid_length, sizeof(int));
  num_grid=(int*) calloc(1, sizeof(int));
  ProbSampleReplace(grid_length, a, perm, 1, num_grid);
  alpha=grid0[num_grid[0]-1];
  alpha=1/alpha-1;
  free(a);
  free(perm);
  free(num_grid);
  free(grid0);
  return alpha;
}

/**************************************************************************/

double nuj_gdp_update(double sigmaj, double alphaj,double *betaj, int j, int m, int grid_length)
{
  double *em, *grid0;
  int *perm, *num_grid;
  double prod_beta, sume, nu;
  int i, k, s;
  grid0=(double *) calloc(grid_length, sizeof(double));
  for (i=0;i<grid_length; i++) grid0[i]=0.00001+i*(1-0.00002)/grid_length;
  /*print_matrix("grid0=", grid_length, 1, grid0,1);*/
  
  prod_beta=1.0;
  sume=0;
  em=(double *) calloc(grid_length, sizeof(double));
  for (i=0; i<grid_length; i++)
  {
   for (k=0; k<j; k++) 
   {
     for (s=1; s<m; s++) prod_beta=prod_beta*pow( 1+fabs(betaj[k*m+s]-betaj[k*m+s-1])*grid0[i]/(sqrt(sigmaj)*(1-grid0[i])),-(alphaj+1)) * grid0[i]/(1-grid0[i]);
   }
   em[i]=prod_beta;
   /*Rprintf("prod_beta=%1.10f\n", prod_beta);*/
   prod_beta=1;
   if(!R_FINITE(em[i])||ISNAN(em[i])||ISNA(em[i])) em[i]=0;
   sume=sume+em[i];
  }
  /*Rprintf("sume=%1.10f\n", sume);*/
  
  if (sume>0) {for (i=0; i<grid_length; i++) em[i]=em[i]/sume;} else {for (i=0; i<grid_length; i++) em[i]=1/grid_length;}
  
  /*print_matrix("em from fct nuj_gdp_update", grid_length,1,em,1 );*/
  perm=(int*) calloc(grid_length, sizeof(int));
  num_grid=(int*) calloc(1, sizeof(int));
  ProbSampleReplace(grid_length, em, perm, 1, num_grid);
  nu=grid0[num_grid[0]-1];
  nu=1/nu-1;
  free(em);
  free(perm);
  free(num_grid);
  return nu;
}


/**************************************************************************/
double lambdakjs_gdp_update(double alphaj, double nuj, double sigmaj, double betakjs, double betakjs1)
{
  double scl, shp, lambdakjs;
  scl=fabs(betakjs-betakjs1)/sqrt(sigmaj)+nuj;
  scl=1/scl;/**/
  shp=alphaj+1;
  lambdakjs=rgamma(shp, scl);
  return lambdakjs; 
}
/**************************************************************************/
double taukjs_gdp_update(double sigmaj, double lambdakjs, double betakjs, double betakjs1)
{
  double taukjs, mu, rho, y2, z, u;
  
  mu=sqrt(sigmaj)*fabs(lambdakjs/(betakjs-betakjs1));
  rho=pow(lambdakjs, 2);
  y2=rchisq(1);
  z=mu+(pow(mu,2)*y2-mu*sqrt(4*mu*rho*y2+pow(mu*y2, 2)))/(2*rho);
  u=runif(0,1);
  if (u<=mu/(mu+z)) taukjs=1/z; else taukjs=z/pow(mu, 2);
  return taukjs; 
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
void post_betakj_gdp(double *Bkj_mat, double *taukj, double *yj, double *muj_k, double sigj, int m, int n, double *meankj,  double *half_covkj)
/*
* this function computes the peoterior mean and covariance for beta_j
*  meanj is a mj vector and covj,half_covj mj x mj matrices
* 
*/
{
  double *BtB, *Bt, *DtT1D, *T1diagcov, *yj0, *covkj, *Dmat;
  int i,get_inv, m1;
  m1=m-1;
  BtB=(double *) calloc(m*m, sizeof(double));
  Bt=(double *) calloc(m*n, sizeof(double));
  yj0=(double *) malloc(n*sizeof(double));
  covkj=(double*)calloc(m*m, sizeof(double));
  DtT1D=(double *) calloc(m*m, sizeof(double));
  Dmat=(double *) calloc((m-1)*m, sizeof(double));
  for (i=0; i<m-1; i++) {Dmat[(m-1)*i+i]=-1; Dmat[(m-1)*(i+1)+i]=1;}
  /*print_matrix("Dmat", m-1,m, Dmat, m-1);*/

  T1diagcov=(double *) calloc(m-1, sizeof(double));
  for (i=0; i<(m-1); i++) T1diagcov[i]=1/taukj[i];
  /*print_matrix("taukj=", m-1,1,taukj,1);*/
  
  getXtWX(DtT1D, Dmat,T1diagcov,&m1,&m,covkj);
  /*print_matrix("DtT1D", m,m, DtT1D, m);*/
  memset(covkj, 0, m*m*sizeof(double));
  free(Dmat);
 
  getXtX(Bkj_mat,&n,&m, BtB);
  for (i=0; i<m*m; i++) BtB[i]=(BtB[i]+DtT1D[i]);
  /*print_matrix("BtB", m,m, BtB, m);
  Rprintf("\n");*/
  get_inv=1;
  T1diagcov=(double *) realloc(T1diagcov, m*sizeof(double));
  qr_ldet_inv(BtB,&m,covkj,&get_inv, T1diagcov);/**/
   for (i=0; i<m*m; i++) half_covkj[i]=covkj[i]; /*covkj is now computed and half_covj is currently equal to covj*/
  free(T1diagcov);
  free(DtT1D);
  
  BtB=(double *) realloc(BtB, m*n*sizeof(double));
  transpose_mat(Bkj_mat, m, n, Bt);  /*Bt now contains t(Bj_mat)*/
  matprod(covkj, m, m, Bt,  n, m, BtB); /*BtB now contains (covkj)^{-1}%*%t(Bj_mat)*/
  for (i=0; i<n; i++) yj0[i]=yj[i]-muj_k[i];
  matprod(BtB,  n, m, yj0, 1, n, meankj); /*meanj now contains (BtB)^{-1}%*%t(Bj_mat)%*%yj */
  /*print_matrix("meanj", m,1, meankj, 1);*/
  free(Bt);
  
  mroot(half_covkj,&m,&m);/*now half_covj contains (covj)^0.5 but stored by rows*/
  BtB=(double *) realloc(BtB, m*m*sizeof(double));
  transpose_mat(half_covkj, m, m, BtB);      /*now BtB contains (covj)^0.5 stored in a right way: by columns*/
  for (i=0; i<m*m; i++) half_covkj[i]=BtB[i]*sqrt(sigj);/*now half_covj contains (covj)^0.5 stored in a right way: by columns*/
  free(BtB);
  free(yj0);
  free(covkj);
}

/*****************************************************************************************************/

void one_step_gibbs_gdp(int grid_length, int m, int n, int p, double *Ball_mat, double *data, double *mu, double *all_beta, double *sigma_vec, double *tau, double *alpha, double *nu)
{
 int j,i,s, mj1, k;
 double sig_j, taukjs, nuj, alphaj, lambdakjs;
 double *betaj,*betakj, *half_covkj, *yj, *muj, *muj_k, *mu_betakj, *tauj, *taukj, *Bkj_mat;

 yj=(double*)malloc(n*sizeof(double));
 muj=(double*)malloc(n*sizeof(double)); 
 muj_k=(double*)malloc(n*sizeof(double));
 half_covkj=(double*)malloc(m*m*sizeof(double));
 mu_betakj=(double*)malloc(m*sizeof(double));
 betakj=(double*)malloc(m*sizeof(double));
 taukj=(double*)malloc((m-1)*sizeof(double)); 
 Bkj_mat=(double *) calloc(n*m,sizeof(double));
 /*print_matrix("begining of one step tau=",p*(p+1)*(m-1)/2,1 ,tau, 1);*/
 
 for (j=0; j<p; j++)
  { 
    mj1=m*(j+1);
    for (i=0; i<n; i++) yj[i]=data[n*j+i];
    for (i=0; i<n; i++) muj[i]=mu[j*n+i];
    betaj=(double*)malloc(mj1*sizeof(double));
    for (s=0; s<mj1; s++) betaj[s]=all_beta[j*(j+1)*m/2+s];
    tauj=(double*)malloc((m-1)*(j+1)*sizeof(double));
    for (s=0; s<(m-1)*(j+1); s++) tauj[s]=tau[j*(j+1)*(m-1)/2+s];
    /*print_matrix("tauj transmitted to sigmaj_gdp_update:",(m-1)*(j+1),1, tauj, 1 );*/
    sig_j=sigmaj_gdp_update(yj, muj, n, j, m, tauj, betaj);
    Rprintf("sig_j=%f\n", sig_j);/**/
    /*sig_j=10;*/
    sigma_vec[j]=sig_j;
    /*nuj=nu[j];
    alphaj=alphaj_gdp_update(sig_j, nuj,betaj, j+1, m, grid_length);
    alpha[j]=alphaj;
    nuj=nuj_gdp_update(sig_j, alphaj, betaj, j+1, m, grid_length);
    nu[j]=nuj;
    Rprintf("alphaj=%f,  nuj=%f\n", alphaj, nuj);*/
    nuj=1; alphaj=1;   /**/
    for (k=0; k<j+1; k++)
    {
      for (s=0; s<m; s++) betakj[s]=betaj[k*m+s];
      for (s=0; s<m-1; s++) 
      {
	lambdakjs=lambdakjs_gdp_update(alphaj, nuj, sig_j, betakj[s], betakj[s+1]);
	/*Rprintf("lambdakjs=%f\n",lambdakjs);*/
	taukjs=taukjs_gdp_update(sig_j, lambdakjs, betakj[s], betakj[s+1]);
	/*Rprintf("taukjs=%f\n",taukjs);*/
	taukj[s]=taukjs; 
	tau[j*(j+1)*(m-1)/2+k*(m-1)+s]=taukjs;
      }  
      for (s=0; s<m*n; s++) Bkj_mat[s]=Ball_mat[m*n*k+s];
      memset(muj_k, 0, n*sizeof(double));
      matprod(Bkj_mat, m, n, betakj, 1, m, muj_k);
      for (s=0; s<n; s++) muj_k[s]=muj[s]-muj_k[s];   
      post_betakj_gdp(Bkj_mat, taukj, yj, muj_k, sig_j, m, n, mu_betakj, half_covkj);
      sim_betakj_psplines(m, mu_betakj, half_covkj, betakj);
      for (s=0; s<m; s++) all_beta[j*(j+1)*m/2+k*m+s]=betakj[s];
    }
    free(betaj);
    free(tauj);
  }
  /*print_matrix("end of one step  tau=",p*(p+1)*(m-1)/2,1 ,tau, 1);*/
  
 free(yj);
 free(muj);
 free(betakj);
 free(muj_k);
 free(half_covkj);
 free(mu_betakj);
 free(taukj);
 free(Bkj_mat);
 }
/*****************************************************************************************************/

SEXP gibbs_BALVM_GDP(SEXP Rtot_draws, SEXP Rknots, SEXP RNknots, SEXP RknotsI, SEXP Rorder,
          SEXP Rn, SEXP Rm, SEXP Rp, SEXP RGrid_length, SEXP Rgrid, SEXP Rdata,
		    SEXP Reta_start, SEXP Rallbeta_start, SEXP Rdens_comp, SEXP Rdrws_step)
{ 
 int i, j, m, p, n, order, grid_length;
 int pn, t, tot_draws, Nknots;
 int dens_comp, drws_step, it_step, eff_draws;
 
 double *data, *ystar, *knots, *knotsI;
 double *eta, *all_beta;
 double *mu, *Ball_mat;
 double  *tau, *alpha, *nu;
 double *sigma_vec;
 double *grid; 
 double ***Pijs;

 /* VARIABLES FOR COMPUTING DENSITY ESTIMATION */
 double *eta_spec, *allbeta_spec, *Ball_mat_spec, *mu_spec;
 
 /* FILE NAMES FOR WRITING RESULTS*/
 char f_mu_name[13]="mu_draws_gdp", f_sig_name[16]="sigma_draws_gdp", 
      f_eta_name[14]="eta_draws_gdp", f_mu_spec_name[18]="mu_spec_draws_gdp", 
      f_sig_spec_name[21]="sigma_spec_draws_gdp",f_beta_name[16]="beta_draws_gdp", f_allmu_name[17]="allmu_draws_gdp";
      
 FILE *out_mu, *out_sigma, *out_eta, *out_beta,
      *out_mu_spec, *out_sigma_spec,
      *in_mu, *in_sigma, *in_eta,*in_beta, 
      *in_mu_spec, *in_sigma_spec, *out_allmu, *out_allmu_spec;
 
 SEXP ans, Rmu, Rsigma, Rsigma_spec, Reta, Rbeta;

 Nknots=INTEGER(RNknots)[0];
 drws_step=INTEGER(Rdrws_step)[0];
 m=INTEGER(Rm)[0];
 p=INTEGER(Rp)[0];
 order=INTEGER(Rorder)[0];
 n=INTEGER(Rn)[0];
 pn=p*n;
 grid_length=INTEGER(RGrid_length)[0];
 tot_draws=INTEGER(Rtot_draws)[0];
 dens_comp=INTEGER(Rdens_comp)[0];

 Pijs=DM3D(n,p,grid_length);
 
 data=(double *) malloc(pn*sizeof(double));
 eta=(double *) malloc(pn*sizeof(double));
 
 for (i=0; i<pn; i++) {data[i]=REAL(Rdata)[i]; eta[i]=REAL(Reta_start)[i];}

 grid=(double *) malloc(grid_length*sizeof(double));
 for (i=0; i<grid_length; i++) grid[i]=REAL(Rgrid)[i];

 all_beta=(double *) malloc(p*(p+1)*m*sizeof(double)/2);
 for (i=0; i<p*(p+1)*m/2; i++) all_beta[i]=REAL(Rallbeta_start)[i];
 print_matrix("all_beta=", p*(p+1)*m/2, 1,all_beta, p*(p+1)*m/2);

 knots=(double *) malloc(Nknots*sizeof(double));
 for (i=0; i<Nknots; i++) knots[i]=REAL(Rknots)[i];
 knotsI=(double *) malloc((Nknots-1)*sizeof(double));
 for (i=0; i<(Nknots-1); i++) knotsI[i]=REAL(RknotsI)[i];

 mu=(double *) malloc(pn*sizeof(double));
 Ball_mat=(double *) calloc(n*m*p,sizeof(double));

 out_allmu=fopen(f_allmu_name, "wb");
 mu_B_fct (n, p, order, m, Nknots, knots, knotsI,  all_beta, eta, Ball_mat, mu,out_allmu);
 /*print_matrix("mu=", n, p,mu, n);
 print_matrix("Ball_mat=", n, p*m,mu, n);*/

 sigma_vec=(double *) malloc(p*sizeof(double));
 ystar=(double *) malloc(grid_length*p*sizeof(double));
 tau=(double*)malloc(p*(p+1)*(m-1)*sizeof(double)/2);
 for (i=0; i<p*(p+1)*(m-1)/2; i++) tau[i]=1; 
 alpha=(double*)malloc(p*sizeof(double));
 nu=(double*)malloc(p*sizeof(double));
 for (i=0; i<p; i++) {alpha[i]=1; nu[i]=1;} /* initial values for nu and alpha*/


  if (dens_comp==1)
  {
   mu_spec=(double *) calloc(p*grid_length, sizeof(double));
   Ball_mat_spec=(double *) calloc(grid_length*m*p,sizeof(double));
   eta_spec=(double *) malloc(p*grid_length*sizeof(double));
   for (j=0; j<p; j++) {for (i=0; i<grid_length; i++) eta_spec[j*grid_length+i]=grid[i];}
   allbeta_spec=(double *) malloc(p*(p+1)*m*sizeof(double)/2);
   for (i=0; i<p*(p+1)*m/2; i++) allbeta_spec[i]=all_beta[i];
   mu_B_fct (grid_length, p, order, m, Nknots, knots, knotsI,  allbeta_spec, eta_spec, Ball_mat_spec, mu_spec,out_allmu_spec);
   out_mu_spec=fopen(f_mu_spec_name, "wb");
   out_sigma_spec=fopen(f_sig_spec_name, "wb");
  }

 out_mu=fopen(f_mu_name, "wb");
 out_sigma=fopen(f_sig_name, "wb");
 out_eta=fopen(f_eta_name, "wb");
 out_beta=fopen(f_beta_name, "wb");

 GetRNGstate();
 
 /*starting draws*/
 
 it_step=0;
 
for (eff_draws=0; eff_draws<tot_draws; )
 {
 one_step_gibbs_gdp(grid_length, m, n, p, Ball_mat, data, mu, all_beta, sigma_vec, tau, alpha, nu);
 comp_etaY_grid(n, p, grid_length, order, m, Nknots, knots, knotsI, grid, all_beta, sigma_vec, data, eta, ystar, Pijs, dens_comp);
 mu_B_fct (n, p, order, m, Nknots, knots, knotsI,  all_beta, eta, Ball_mat, mu,out_allmu); 
 it_step=it_step+1;
 
 if (it_step==drws_step)
 {
  fsafewrite(mu, sizeof(double), pn, out_mu);
  fsafewrite(sigma_vec, sizeof(double),p, out_sigma);
  fsafewrite(eta, sizeof(double), pn, out_eta);
  fsafewrite(all_beta, sizeof(double), p*(p+1)*m/2, out_beta);
  if (dens_comp==0)
   {
     it_step=0;
     eff_draws=eff_draws+1;

   }
 }

Rprintf("it_step=%d,  eff_draws=%d\n",it_step, eff_draws);
}


 PutRNGstate();
 fclose(out_allmu);
 fclose(out_mu);
 fclose(out_sigma);
 fclose(out_eta);
 fclose(out_beta);

 if(dens_comp==1)
 {
  fclose(out_mu_spec);
  fclose(out_sigma_spec);
 }/**/

 free(data);
 free(grid);
 free(knots);
 free(knotsI);
 free(Ball_mat);
 free(ystar);
 Dfree3d(Pijs);
 free(tau);
 free(alpha);
 free(nu);

 if (dens_comp==1)
  {
  free(eta_spec);
  free(allbeta_spec);
  free(Ball_mat_spec);
  }
 
 in_mu=fopen(f_mu_name, "rb");
 in_sigma=fopen(f_sig_name, "rb");
 in_eta=fopen(f_eta_name, "rb");
 in_beta=fopen(f_beta_name, "rb");
 
 if (dens_comp==1)
  {
    in_mu_spec=fopen(f_mu_spec_name, "rb");
    in_sigma_spec=fopen(f_sig_spec_name, "rb");
  }

  if (dens_comp==1)
  {
   PROTECT(ans = allocVector(VECSXP, 3*tot_draws+2));
   PROTECT(Rsigma = allocMatrix(REALSXP,  p, tot_draws));
   PROTECT(Rsigma_spec = allocMatrix(REALSXP,  p, tot_draws));
  } else
    
  {
   PROTECT(ans = allocVector(VECSXP, 3*tot_draws+1));
   PROTECT(Rsigma = allocMatrix(REALSXP,  p, tot_draws));
  }
 for (t=0; t<tot_draws; t++) {
    fsaferead(mu, sizeof(double), pn, in_mu);
    fsaferead(sigma_vec, sizeof(double), p, in_sigma);
    fsaferead(eta, sizeof(double), pn, in_eta);
    fsaferead(all_beta, sizeof(double), p*(p+1)*m/2, in_beta);
    PROTECT(Rmu = allocMatrix(REALSXP, n, p));
    PROTECT(Reta = allocMatrix(REALSXP, n, p));
    PROTECT(Rbeta = allocVector(REALSXP, p*(p+1)*m/2));
    for(i = 0; i < pn; i++) {REAL(Rmu)[i] = mu[i]; REAL(Reta)[i] = eta[i];}
    for(i = 0; i < p*(p+1)*m/2; i++) {REAL(Rbeta)[i] = all_beta[i];}
    for(i = 0; i < p; i++) REAL(Rsigma)[t*p+i] = sigma_vec[i];
    SET_VECTOR_ELT(ans, t, Rmu);
    SET_VECTOR_ELT(ans, tot_draws+t, Reta);
    SET_VECTOR_ELT(ans, 2*tot_draws+t, Rbeta);
    UNPROTECT(3);
    
    if (dens_comp==1)
    {
      fsaferead(mu_spec, sizeof(double), p*grid_length, in_mu_spec);
      fsaferead(sigma_vec, sizeof(double), p, in_sigma_spec);
      PROTECT(Rmu = allocMatrix(REALSXP, grid_length, p));
      for(i = 0; i < p*grid_length; i++) REAL(Rmu)[i] = mu_spec[i];
      for(i = 0; i < p; i++) REAL(Rsigma_spec)[t*p+i] = sigma_vec[i];
      SET_VECTOR_ELT(ans, 2*tot_draws+t, Rmu);
      UNPROTECT(1);
    }
 }
 if (dens_comp==1) {
  SET_VECTOR_ELT(ans, 3*tot_draws, Rsigma);
  SET_VECTOR_ELT(ans, 3*tot_draws+1, Rsigma_spec);
  UNPROTECT(3);
  fclose(in_mu);
  fclose(in_sigma);
  fclose(in_mu_spec);
  fclose(in_sigma_spec);
  free(mu);
  free(mu_spec);
  free(sigma_vec);
 } else {
   SET_VECTOR_ELT(ans, 3*tot_draws, Rsigma);
   UNPROTECT(2);
    fclose(in_mu);
    fclose(in_sigma);
    fclose(in_eta);
    fclose(in_beta);
    free(mu);
    free(sigma_vec);
    free(eta);
    free(all_beta);
 }
 return (ans);
}

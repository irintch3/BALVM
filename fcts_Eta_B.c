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
    Rprintf( " %s", desc );
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
    if(n <= SMALL) {
	/* might do this repeatedly, so speed matters */
	HL = (int *)alloca(n * sizeof(int));
	q = (double *) alloca(n * sizeof(double));
	R_CheckStack();
    } else {
	/* Slow enough anyway not to risk overflow */
	HL = Calloc(n, int);
	q = Calloc(n, double);
    }
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

/************************************************************************/

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
/*
This procedure computes and simulates the new eta and ystar
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


SEXP splb_check(SEXP Rknots, SEXP Rnknots, SEXP RknotsI, SEXP Rorder, SEXP Reta, SEXP Rm,
                SEXP Rallbeta, SEXP Rn, SEXP Rp, SEXP RS, SEXP RG, SEXP Rall_sigma, SEXP Rdata)
{
/* evaluate the non-zero B-spline basis functions (or their derivatives)
 * at xvals.  */
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


void mu_B_fct (int n, int p, int order, int m, int nknots, double *knots, double *knotsI,  double *all_beta, double *eta, double *Ball_mat, double *mu)
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


double sim_sigma_j(double *yj, double *muj, int n, double asigj, double bsigj)
{
  double sigj, sclsigj, shpsigj;
  int i;
  sclsigj=1/bsigj;
  for (i=0; i<n; i++) sclsigj=sclsigj+0.5*pow(yj[i]-muj[i],2);
  shpsigj=asigj+0.5*n;
  sigj=1/rgamma(shpsigj, sclsigj);
  return sigj;
}

double sim_taukj_naive(double ataukj, double btaukj, double *betakj, int m)
{
  double taukj, scltaukj, shptaukj;
  int i;
  scltaukj=1/btaukj;
  for (i=0; i<m; i++) scltaukj=scltaukj+0.5*pow(betakj[i],2);
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

void post_betaj_mean_cov(double *Bj_mat, double *tauj, double *yj, double sigj, int j, int m, int n, double *meanj, double *covj, double *half_covj) 
/*
* this function computes 
*
*/
{
  double *BtB, *Bt, *sigtauj1, *diagcovj;
  int i,k, mj1, get_inv;
  mj1=m*(j+1);
  BtB=(double *) calloc(mj1*mj1, sizeof(double));
  Bt=(double *) calloc(mj1*n, sizeof(double));
  sigtauj1=(double *) calloc((j+1), sizeof(double));
  diagcovj=(double *) calloc(mj1, sizeof(double));
  getXtX(Bj_mat,&n,&mj1, BtB);
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
  free(BtB);/*don't need BtB anymore will use it for other stuff*/ 
  free(sigtauj1);
  free(diagcovj);
  
  BtB=(double *) calloc(mj1*n, sizeof(double));
  transpose_mat(Bj_mat, mj1, n, Bt); /*Bt now contains t(Bj_mat)*/
  matprod(covj, mj1, mj1, Bt,  n, mj1, BtB); /*BtB now contains (BtB)^{-1}%*%t(Bj_mat)*/
  matprod(BtB,  n, mj1, yj, 1, n, meanj); /*meanj now contains (BtB)^{-1}%*%t(Bj_mat)%*%yj */
  /*print_matrix("meanj", n,1, meanj, 1);*/
  free(BtB); /*don't need BtB anymore*/
  free(Bt);
  
  mroot(half_covj,&mj1,&mj1);/*now half_covj contains (covj)^0.5 but stored by rows*/
  BtB=(double *) calloc(mj1*mj1, sizeof(double));
  transpose_mat(half_covj, mj1, mj1, BtB);
  for (i=0; i<mj1*mj1; i++) half_covj[i]=BtB[i];/**/
  free(BtB);
}

/****************************************************************************************************************************************************/

void sim_betaj(int j, int m, double *meanj, double *half_covj, double *betaj)
{
  int k, mj1;
  double *hlp;
  mj1=m*(j+1);
  hlp=(double *) calloc(mj1, sizeof(double));
  for (k=0; k<mj1; k++) betaj[k]=rnorm(0, 1);
  matprod(betaj, 1, mj1, half_covj, mj1, mj1, hlp);
  for (k=0; k<mj1; k++) betaj[k]=hlp[k]+meanj[k];
  free(hlp);
}

/****************************************************************************************************************************************************/

SEXP betaj_post_check(SEXP Rknots, SEXP Rnknots, SEXP RknotsI, SEXP Rorder, SEXP Reta, SEXP Rm,
                SEXP Rallbeta, SEXP Rn, SEXP Rp, SEXP RS, SEXP RG, SEXP Rall_sigma, SEXP Rdata, SEXP Rtau)
{
/* this function checks whether we compute correctly parameters for betaj posterior distribution*/
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

 comp_etaY_grid (n, p, grid_length, order, m, Nknots, knots, knotsI, grid, all_beta, sigma_vec, data, eta, ystar, Pijs);/**/

} /* end draws*/


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

SEXP B_mu_check(SEXP Rknots, SEXP Rnknots, SEXP RknotsI, SEXP Rorder, SEXP Reta, SEXP Rm,
                SEXP Rallbeta, SEXP Rn, SEXP Rp, SEXP RS, SEXP RG, SEXP Rall_sigma, SEXP Rdata)
{
/* evaluate the non-zero B-spline basis functions (or their derivatives)
 * at xvals.  */
    int nknots, order, m, n, i, p, nbeta, S/*, s,s1,j*/ ;
    double *knots, *knotsI, *eta, *all_beta;
    SEXP ans, RBall_mat, Rmu;/**/
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

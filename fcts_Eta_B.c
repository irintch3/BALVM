/*  Routines for manipulating B-splines are by Copyright (C) 1998 Douglas M. Bates and William N. Venables.
 *
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


/*splines manipulating fct 
* set sp->curs to the index of the first knot position > x.
* Special handling for x == sp->knots[sp->nknots - sp-order + 1] */

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


/*splines manipulating fct 
* set sp->rdel and sp->ldel 
* as fct of point x and 
*/
static void
diff_table(splPTR *sp, double x, int ndiff)
{
  int i;
  for (i = 0; i < ndiff; i++) {
      sp->rdel[i] = sp->knots[sp->curs + i] - x;
      sp->ldel[i] = x - sp->knots[sp->curs - (i + 1)];
  }
}

/* fast evaluation of spline basis functions */
static void basis_funcs(splPTR *sp, double x, double *b)
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


/* 
 * evaluates the non-zero B-spline basis functions (and their derivatives)
 * at xvals (of length nx) with *knots of length nk with derivs (of length nd)
 * the result is in *val
 */
 
void spline_basis(double *knots, int nk, int order, double *xvals, int nx, int *derivs, int nd, int *offsets, 
                  double *val)
{
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

/*
* print matrix contained in *a having m columns and n rows
* with description in desc
*/
void print_matrix( char* desc, int m, int n, double* a, int lda )
{
    int i, j;
    Rprintf( " %s", desc );
    for( i = 0; i < m; i++ ) {
            for( j = 0; j < n; j++ ) Rprintf( " %6.9f", a[i+j*lda] );
            Rprintf( "\n" );
    }
}

/*
* selects the subvector od *vec1 starting by position a and ending with position b
* without checking whether the positions are smaller than length of *vec1
* resulting vector is in *resvec
*/
void subvec(double *vec1, int a, int  b, double *resvec)
{
  int i=0, length2=b-a;
  for(i = 0; i < length2; i++) resvec[i]=vec1[a+i];
}

/*
 * safely writes *ptr of length n into a FILE 
 */

long fsafewrite(ptr,size,n,stream)
     double *ptr;size_t size;long n;FILE *stream;

{ long i,j,k=0L;
  for (i=0;i<(n/32000L);i++)
  k+=fwrite(ptr+i*32000L,size,(size_t)32000L,stream);
  j=n%32000L;
  k+=fwrite(ptr+i*32000L,size,(size_t)j,stream);
  return(k);
}

/*
 * safely reads *ptr of length n from a FILE 
 */
long fsaferead(ptr,size,n,stream)
     double *ptr;size_t size;long n;FILE *stream;

{ long i,j=32000L,k=0L;
  for (i=0;i<(n/32000L);i++)
  k+=fread(ptr+i*32000L,size,(size_t)j,stream);
  j=n%32000L;
  k+=fread(ptr+i*32000L,size,(size_t)j,stream);
  return(k);
}

/*
 * transposes *mat  
 */
void transpose_mat(double *mat, int ncol, int nrow, double *resmat)
{
 int j=0, i=0;
       for(i = 0; i < nrow; i++)    /* column number for the element of resmat */
         for (j = 0; j < ncol; j++)     /* line number for the element of resmat */
           {
              resmat[i* ncol+j]=mat[j* nrow+i];
           }
}

/*
 * computes the procuct of mat1 by mat2 
 */
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

/* A  version unequal probability sampling using Walker's alias method, based on Alg 3.13B in
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

  /*
  computes multivarite normal density
  all_sigma is vector containing the diagonal elements of the covariance matrix
  y is the point where the density is computed
  u is the location parameter
  */

double dmvnorm(int dim, double *all_sigma, double *y, double *mu)
{
  int i;
  double val;
  double sum_sq, lndet;
  sum_sq=0;
  lndet=0;
  
  for (i=0; i<dim; i++){sum_sq=sum_sq+pow(y[i]-mu[i], 2.0)/(2*all_sigma[i]); lndet=lndet+0.5*log(all_sigma[i]*2*PI);}
  val=exp(-sum_sq-lndet);
  return val;
}


/*
 * allocates memory to a three dimensional array
 */
 
double ***DM3D(int n,int row,int col) 
{
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
/*
 * frees memory taken by a three dimensional array 
 */
void Dfree3d(double ***pa) 
{

/* free the data */
free(**pa);
/* free the row pointers */
free(*pa);
/* free the 3D pointers */
free(pa);
}
/*
 * allocates memory to a two dimensional array
 */
double **DM2D(int n,int row) 
{
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
/*
 * frees memory taken by a two dimensional array 
 */
void Dfree2d(double **pa) 
{

free(*pa);
free(pa);
}


/*
 * samples sigma_j 
 */
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

/*
* samples tau_kj
*/
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


/*
 * this function computes meanj and halfcovj necessary for update of betaj
 */

void post_betaj_mean_cov(double *Bj_mat, double *tauj, double *yj, double sigj, int j, int m, int n, 
                         double *meanj, double *covj, double *half_covj) 
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
  
  for (i=0; i<(j+1); i++) 
  {
    for (k=0; k<m;k++) BtB[mj1*m*i+mj1*k+m*i+k]=BtB[mj1*m*i+mj1*k+m*i+k]+sigtauj1[i];
  }
  
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
  
  mroot(half_covj,&mj1,&mj1);
  /*now half_covj contains (covj)^0.5 but stored by rows*/
  BtB=(double *) calloc(mj1*mj1, sizeof(double));
  transpose_mat(half_covj, mj1, mj1, BtB);
  for (i=0; i<mj1*mj1; i++) half_covj[i]=BtB[i];
  free(BtB);
}


/*
* sample betaj
*/
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


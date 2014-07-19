
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
    int order,  		/* order of the spline */
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

static int set_cursor(splPTR *sp, double x);


/*splines manipulating fct 
* set sp->rdel and sp->ldel 
* as fct of point x and 
*/
static void
diff_table(splPTR *sp, double x, int ndiff);

/* fast evaluation of spline basis functions */

static void 
basis_funcs(splPTR *sp, double x, double *b);


/* "slow" evaluation of (derivative of) basis functions */
static double evaluate(splPTR *sp, double x, int nder);



/* 
 * evaluates the non-zero B-spline basis functions (and their derivatives)
 * at xvals (of length nx) with *knots of length nk with derivs (of length nd)
 * the result is in *val
 */
 
void spline_basis(double *knots, int nk, int order, double *xvals, int nx, int *derivs, int nd, int *offsets, 
                  double *val);

/*
* print matrix contained in *a having m columns and n rows
* with description in desc
*/
void print_matrix( char* desc, int m, int n, double* a, int lda );

/*
* selects the subvector od *vec1 starting by position a and ending with position b
* without checking whether the positions are smaller than length of *vec1
* resulting vector is in *resvec
*/
void subvec(double *vec1, int a, int  b, double *resvec);


/*
 * safely writes *ptr of length n into a FILE 
 */

long fsafewrite(double *ptr, size_t size, long n, FILE *stream);


/*
 * safely reads *ptr of length n from a FILE 
 */
long fsaferead(double *ptr, size_t size, long n, FILE *stream);


/*
 * transposes *mat  
 */
void transpose_mat(double *mat, int ncol, int nrow, double *resmat);

/*
 * computes the procuct of mat1 by mat2 
 */
void matprod(double *mat1, int ncol1, int nrow1, double *mat2, int ncol2, int nrow2, double *resmat);


/*
 *  Unequal Probability Sampling.
 *
 *  Modelled after Fortran code provided by:
 *    E. S. Venkatraman <venkat@biosta.mskcc.org>
 *  but with significant modifications in the
 *  "with replacement" case.
 */

/* Unequal probability sampling; with-replacement case */

void ProbSampleReplace(int n, double *p, int *perm, int nans, int *ans);


/* A  version unequal probability sampling using Walker's alias method, based on Alg 3.13B in
   Ripley (1987).
 */
#define SMALL 10000
static void walker_ProbSampleReplace(int n, double *p, int *a, int nans, int *ans);

  /*
  computes multivarite normal density
  all_sigma is vector containing the diagonal elements of the covariance matrix
  y is the point where the density is computed
  u is the location parameter
  */

double dmvnorm(int dim, double *all_sigma, double *y, double *mu);


/*
 * allocates memory to a three dimensional array
 */
 
double ***DM3D(int n,int row,int col);

/*
 * frees memory taken by a three dimensional array 
 */
void Dfree3d(double ***pa);

/*
 * allocates memory to a two dimensional array
 */
double **DM2D(int n,int row); 

/*
 * frees memory taken by a two dimensional array 
 */
void Dfree2d(double **pa); 



/*
 * samples sigma_j 
 */
double sim_sigma_j(double *yj, double *muj, int n, double asigj, double bsigj);


/*
* samples tau_kj
*/
double sim_taukj_naive(double ataukj, double btaukj, double *betakj, int m);



/*
 * this function computes meanj and halfcovj necessary for update of betaj
 */

void post_betaj_mean_cov(double *Bj_mat, double *tauj, double *yj, double sigj, int j, int m, int n, 
                         double *meanj, double *covj, double *half_covj); 


/*
* sample betaj
*/
void sim_betaj(int j, int m, double *meanj, double *half_covj, double *betaj);



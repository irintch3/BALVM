#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#define LOW -1.0e15

/* call LA_PACK to get pivoted QR decomposition of x
   tau is an array of length min(r,c)
   pivot is array of length c, zeroed on entry, pivoting order on return.
   On exit upper triangle of x is R. Below upper triangle plus tau
   represent reflectors making up Q.
   pivoting is always performed (not just when matrix is rank deficient), so
   leading diagonal of R is in descending order of magnitude.
   library(mgcv)
   r<-4;c<-3
   X<-matrix(rnorm(r*c),r,c)
   pivot<-rep(1,c);tau<-rep(0,c)
   um<-.C("mgcv_qr",as.double(X),as.integer(r),as.integer(c),as.integer(pivot),as.double(tau))
   qr.R(qr(X));matrix(um[[1]],r,c)[1:c,1:c]
*/

void mgcv_qr(double *x, int *r, int *c,int *pivot,double *tau);


/* applies k reflectors of Q of a QR decomposition to r by c matrix b.
   Apply Q from left if left!=0, right otherwise.
   Transpose Q only if tp!=0.
   Information about Q has been returned from mgcv_qr, and is stored in tau and
   below the leading diagonal of a.
   library(mgcv)
   r<-4;c<-3
   X<-matrix(rnorm(r*c),r,c)
   qrx<-qr(X)
   pivot<-rep(1,c);tau<-rep(0,c)
   um<-.C("mgcv_qr",a=as.double(X),as.integer(r),as.integer(c),as.integer(pivot),tau=as.double(tau))
   y<-1:4;left<-1;tp<-0;cy<-1
   er<-.C("mgcv_qrqy",as.double(y),as.double(um$a),as.double(um$tau),as.integer(r),as.integer(cy),as.integer(c),
        as.integer(left),as.integer(tp),PACKAGE="mgcv")
   er[[1]];qr.qy(qrx,y)

*/
void mgcv_qrqy(double *b,double *a,double *tau,int *r,int *c,int *k,int *left,int *tp);


/* Finds C = R^{-1} B where R is the c by c matrix stored in the upper triangle
   of r by c argument R. B is c by bc. (Possibility of non square argument
   R facilitates use with output from mgcv_qr). This is just a standard back
   substitution loop.
*/
void mgcv_backsolve(double *R,int *r,int *c,double *B,double *C, int *bc);


/* Obtains the log|X| and the inverse of X (r by r), by pivoted QR decomposition.
   The inverse is returned (unpivoted) in Xi.
   The function returns log|X| as its value.
   X is overwirtten in the process
*/
void qr_ldet_inv(double *X,int *r,double *Xi,int *get_inv, double *diagR);



/* forms X'MX as efficiently as possible, where M is a symmetric matrix
   and X is an r by c matrix. X and M are stored column wise.
   work should be an r-vector (longer is no problem).
*/
void getXtMX(double *XtMX,double *X,double *M,int *r,int *c,double *work);





/* a stored in column order, this routine finds the pivoted choleski decomposition of matrix a 
   library(mgcv)
   X<-matrix(rnorm(16),4,4);D<-X%*%t(X)
   rank<-0
   er<-.C("mgcv_chol",as.double(D),as.integer(rep(0,4)),as.integer(4),as.integer(rank))
   rD<-matrix(er[[1]],4,4);piv<-er[[2]]
   chol(D,pivot=TRUE);rD;piv
   n<-length(piv);ind<-1:n;
   ind[piv]<-1:n
   rD<-rD[,ind]
   L<-mroot(D)
   D;t(rD)%*%rD;L%*%t(L)
*/
void mgcv_chol(double *a,int *pivot,int *n,int *rank);


/* finds the minimum rank or supplied rank square root of n by n matrix  A by pivoted choleski 
   decomposition returned in A is B such that B'B = A (this is different from R routine mroot)

   R testing code....
   library(mgcv)
   X<-matrix(rnorm(12),4,3);D<-X%*%t(X)
   rank<-0
   er<-.C("mroot",as.double(D),as.integer(rank),as.integer(4))
   rD<-matrix(er[[1]],er[[2]],4)
   D;t(rD)%*%rD
   
*/
void mroot(double *A,int *rank,int *n);


/* form X'X (nearly) as efficiently as possible 
   r is number of rows, 
   c is number of columns */

void getXtX(double *X,int *r,int *c, double *XtX);


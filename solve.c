#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#define LOW -1.0e15


void mgcv_qr(double *x, int *r, int *c,int *pivot,double *tau)
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
{ int info,lwork=-1,*ip;
  double work1,*work;
  /* workspace query */
  F77_NAME(dgeqp3)(r,c,x,r,pivot,tau,&work1,&lwork,&info);
  lwork=(int)floor(work1);if (work1-lwork>0.5) lwork++;
  work=(double *)calloc((size_t)lwork,sizeof(double));
   /* actual call */
  F77_NAME(dgeqp3)(r,c,x,r,pivot,tau,work,&lwork,&info);
  free(work);
  /*if (*r<*c) lwork= *r; else lwork= *c;*/
  for (ip=pivot;ip < pivot + *c;ip++) (*ip)--;
  /* ... for 'tis C in which we work and not the 'cursed Fortran... */

}


void mgcv_qrqy(double *b,double *a,double *tau,int *r,int *c,int *k,int *left,int *tp)
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

{ char side='L',trans='N';
  int lda,lwork=-1,info;
  double *work,work1;

  if (! *left) { side='R';lda = *c;} else lda= *r;
  if ( *tp) trans='T';
  /* workspace query */
  F77_NAME(dormqr)(&side,&trans,r,c,k,a,&lda,tau,b,r,&work1,&lwork,&info);
  lwork=(int)floor(work1);if (work1-lwork>0.5) lwork++;
  work=(double *)calloc((size_t)lwork,sizeof(double));
  /* actual call */
  F77_NAME(dormqr)(&side,&trans,r,c,k,a,&lda,tau,b,r,work,&lwork,&info);
  free(work);

}

void mgcv_backsolve(double *R,int *r,int *c,double *B,double *C, int *bc)
/* Finds C = R^{-1} B where R is the c by c matrix stored in the upper triangle
   of r by c argument R. B is c by bc. (Possibility of non square argument
   R facilitates use with output from mgcv_qr). This is just a standard back
   substitution loop.
*/
{ int i,j,k;
  double x,*pR,*pC;
  for (j=0;j<*bc;j++) { /* work across columns of B & C */
    for (i = *c-1;i>=0;i--) { /* work up each column of B & C */
      x = 0.0;
      /* for (k=i+1;k<*c;k++) x += R[i + *r * k] * C[k + j * *c]; ...following replaces...*/
      pR = R + i + (i+1) * *r;pC = C + j * *c + i + 1;
      for (k=i+1;k<*c;k++,pR+= *r,pC++) x += *pR * *pC;
      C[i + j * *c] = (B[i + j * *c] - x)/R[i + *r * i];
    }
  }
}


void qr_ldet_inv(double *X,int *r,double *Xi,int *get_inv, double *diagR)
/* Obtains the log|X| and the inverse of X (r by r), by pivoted QR decomposition.
   The inverse is returned (unpivoted) in Xi.
   The function returns log|X| as its value.
   X is overwirtten in the process
*/
{ double *tau,*p,*Qt;
  int *pivot,i,TRUE=1,j;
  /* Allocated working storage ...*/
  pivot = (int *)calloc((size_t)*r,sizeof(int));
  tau = (double *)calloc((size_t)*r,sizeof(double));

  mgcv_qr(X,r,r,pivot,tau); /* get QR=X itself */

  /* evaluate log|X| = sum_i log(|R_ii|) ...*/
  /*/ for (ldet=0.0,p=X,i=0;i<*r;i++,p += *r+1) {ldet += log(fabs(*p));Rprintf("det=%f\n",exp(ldet));};*/
  for (i=0;i<*r;i++) {diagR[i]=fabs(X[i**r+i]); /*/Rprintf("Rii=%f\n",fabs(X[i**r+i]));/*/}
  
  if (*get_inv) {
  /* Now get the inverse of X. X^{-1} = R^{-1}Q' */
    Qt = (double *)calloc((size_t)*r * *r,sizeof(double));
    for (p=Qt,i=0;i<*r;i++,p += *r+1) *p = 1.0;
    mgcv_qrqy(Qt,X,tau,r,r,r,&TRUE,&TRUE); /* Extracting the orthogonal factor Q' */

    mgcv_backsolve(X,r,r,Qt,Xi,r); /* Now Xi contains the row pivoted inverse of X */

    /* Finally unpivot Xi.
       pivot[i] is the unpivoted row that pivoted row i should end up in
    */

    for (p=Xi,j=0;j<*r;j++,p += *r) { /* work column by column using tau as working storage */

      for (i=0;i<*r;i++) tau[pivot[i]] = p[i]; /* ith row of pivoted -> pivot[i] row of unpivoted */
      for (i=0;i<*r;i++) p[i] = tau[i];        /* store unpivoted column in Xi */

    }
    free(Qt);
  } /* end if (*get_inv) */
  free(pivot);free(tau);
  return;
}



void getXtMX(double *XtMX,double *X,double *M,int *r,int *c,double *work)
/* forms X'MX as efficiently as possible, where M is a symmetric matrix
   and X is an r by c matrix. X and M are stored column wise.
   work should be an r-vector (longer is no problem).
*/

{ int i,j;
  double *p,*p1,*p2,*pX0,*pX1,xx,*pM;
  pX0=X;
  for (i=0;i< *c;i++) {
    /* first form MX[:,i] */
    p2 = work + *r;pM=M;
    for (p1=work;p1<p2;pM++,p1++) *p1 = *pX0 * *pM;pX0++;
    for (j=1;j< *r;j++,pX0++)
    for (p1=work;p1<p2;pM++,p1++) *p1 += *pX0 * *pM;
    /* now form ith row and column of X'MX */
    for (pX1=X,j=0;j<=i;j++) {
      for (xx=0.0,p=work;p<p2;p++,pX1++) xx += *p * *pX1;
      XtMX[i * *c + j] = XtMX[j * *c + i] = xx;
    }
  }
}





/* added March, 23, 2012*/

void mgcv_chol(double *a,int *pivot,int *n,int *rank)
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
{ double *work,*p1,*p2,*p;
  int piv=1;
  work=(double *)calloc((size_t) *n,sizeof(double));
  F77_NAME(dchdc)(a,n,n,work,pivot,&piv,rank);
  /* zero stuff below the leading diagonal */
  for (p2=a+ *n,p1=a+1;p2<a+ *n * *n;p1+= *n+1,p2+= *n) for (p=p1;p<p2;p++) *p=0.0;
  free(work);
}

void mroot(double *A,int *rank,int *n)
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
{ int *pivot,erank,i,j;
  double *pi,*pj,*p0,*p1,*B;
  pivot=(int *)calloc((size_t)*n,sizeof(int));
  mgcv_chol(A,pivot,n,&erank);
  if (*rank<=0) *rank = erank;
  /* unscramble the returned choleski factor */ 
  B=calloc((size_t)(*n * *n),sizeof(double));
  /* copy A to B and zero A */
  for (p0=A,p1=B,i=0;i< *n;i++,p0+= *n,p1+= *n)
  for (pi=p0,pj=p1;pi<=p0+i;pi++,pj++) { *pj = *pi; *pi=0.0;}
  /* now copy B into A undoing pivoting */
  for (p0=B,i=0;i< *n;p0+= *n,i++) 
  for (j=pivot[i],pj=p0,pi=A + (j-1)* *n;pj<=p0+i;pi++,pj++) *pi = *pj;
  /* remove trailing rows from choleski factor */
  for (pi=A,p0=A,i=0;i< *n;i++,p0+= *n)
  for (pj=p0;pj<p0+ *rank;pj++,pi++) *pi = *pj; 
  free(pivot);free(B);
}


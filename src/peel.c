#include <R.h>
#include <Rmath.h>
#include <string.h>
#include "type_defs.h"

#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>


/*
#ifndef SINGLE
#endif

#ifdef SINGLE
*/
//#include "mkl_blas.h"
//#include "mkl_lapack.h"
//#endif

////////////////////////////////////////////////////////////////////////////////
//80////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Function declarations for single precision lapack routines that are not 
// declared in R_ext/Lapack.h
////////////////////////////////////////////////////////////////////////////////

/* SGEEV - compute for an N-by-N real nonsymmetric matrix A, the */
/* eigenvalues and, optionally, the left and/or right eigenvectors */
/*
extern void
F77_NAME(sgeev)(const char* jobvl, const char* jobvr,
		const int* n, float* a, const int* lda,
		float* wr, float* wi, float* vl, const int* ldvl,
		float* vr, const int* ldvr,
		float* work, const int* lwork, int* info);
*/
/* SGESV - compute the solution to a real system of linear */
/* equations  A * X = B, */
/*
extern void
F77_NAME(sgesv)(const int* n, const int* nrhs, float* a, const int* lda,
		int* ipiv, float* b, const int* ldb, int* info);
*/

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void init_sm (int n, const double *qmatrix, 
	      double *vals, 
	      double *lv, double *rv, 
	      double *bl, double *br, 
	      double *sd,
	      int *info, int *err)
{
  if (n == 0) return;

  double tol = 1e-15;	// tolerance for zero-check on eigenvalue

  int lwork = -1;
  double amat[n*n];
  double wi[n];

  *info = 0;
  *err = 0;

  memcpy(amat,qmatrix,n*n*sizeof(double));

  {
    double query;

//Call the correct LAPACK routine for the float type requested
#ifdef SINGLE 
    F77_CALL(sgeev)("n","v",&n,amat,&n,vals,wi,
		    lv,&n,rv,&n,
		    &query,&lwork,info);
#else
    F77_CALL(dgeev)("n","v",&n,amat,&n,vals,wi,
		    lv,&n,rv,&n,
		    &query,&lwork,info);
#endif
    
    lwork = (int) query;
  }

  {
    double work[lwork];
//Call the correct LAPACK routine for the float type requested
#ifdef SINGLE
    F77_CALL(sgeev)("n","v",&n,amat,&n,vals,wi,
		    lv,&n,rv,&n,
		    work,&lwork,info);
#else
    F77_CALL(dgeev)("n","v",&n,amat,&n,vals,wi,
		    lv,&n,rv,&n,
		    work,&lwork,info);
#endif
  }

  if (*info != 0) return;

  // find the dominant eigenvalue
  int dom = 0;
  {
    double maxreal = vals[0];
    for (int j = 0; j < n; j++) {
      if (maxreal < vals[j]) {
	maxreal = vals[j];
	dom = j;
      }
    }
  }

  // matrix must be stochastic
  // NB: THIS ASSUMES THERE IS A UNIQUE DOMINANT EIGENVALUE
  // WHEN THE STATIONARY DISTRIBUTION IS NONUNIQUE, THIS IS AN ERROR!
  // FIXME!
  if (vals[dom] > tol) {
     printf("miscreant eigenvalue = %lg\n",vals[dom]);
    *err = 1;
    return;
  }

  // set dominant eigenvalue to zero
  vals[dom] = 0;

  // matrix must have real eigenvalues
  for (int i = 0; i < n; i++) {
    if (fabs(wi[i]) > tol) {
      printf("offending eigenvalue = %lg\n", wi[i]);
      *err = 2;
      return;
    } 
  }

  // normalize stationary distribution
  {
    double sum = 0;
    int i;
    for (i = 0; i < n; i++) sum += rv[i+dom*n];
    for (i = 0; i < n; i++) sd[i] = rv[i+dom*n]/sum;
  }


  // compute left eigenvectors
  {
    int ipvt[n], i, j;

    memcpy(amat,rv,n*n*sizeof(double));
    for (i = 0; i < n*n; i++) lv[i] = 0.0;
    for (i = 0, j = 0; j < n; i += n+1, j++) lv[i] = 1.0; 

//Call the correct LAPACK routine for the float type requested
#ifdef SINGLE
    F77_CALL(sgesv)(&n,&n,amat,&n,ipvt,lv,&n,info);
    if (*info != 0) return;
#else
    F77_CALL(dgesv)(&n,&n,amat,&n,ipvt,lv,&n,info);
    if (*info != 0) return;
#endif
  }

  // compute left and right backward-action matrices
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      bl[i+n*j] = rv[j+n*i];
      br[i+n*j] = lv[j+n*i];
    }
  }

}

void merge_lineages(int n, int nlocus,
		    storage_type *my_data, double *my_charac,
		    const storage_type *her_data, double her_charac) {
  double sum, grandsum;
  storage_type *myp;
  const storage_type *herp;
  int i, j;
  for (j = 0, grandsum = 0, myp = my_data, herp = her_data; j < nlocus; j++) {
    for (i = 0, sum = 0; i < n; i++, herp++) 
      sum += (myp[i] =  myp[i] * *herp);
    if (sum > 0.0) {
      grandsum += log(sum);
      for (i = 0; i < n; i++, myp++) *myp /= sum;
    } else {
      for (i = 0; i < n; i++) myp++; // still advance the pointer in the case of sum == 0
    }
  }
  *my_charac += grandsum + her_charac;
}

void statdist_action (int n, int nlocus, 
		      storage_type *ell, double *charac, 
		      const double *sd, const double *pattern_count) {
  int i, j;
  double sum, grandsum;
  storage_type *ellp;
  for (j = 0, grandsum = 0, ellp = ell; j < nlocus; j++) {
    for (i = 0, sum = 0; i < n; i++) sum += ellp[i]*sd[i];
    grandsum += log(sum*pattern_count[j]);
    for (i = 0; i < n; i++, ellp++) *ellp /= sum;
  }
  *charac += grandsum;
}

void molmodel_action (int n, int nlocus, double t, 
		      storage_type *ell,
		      double relax,
		      const double *vals, 
		      const double *lv, const double *rv) {
  if (t != 0) {
    double phi[n];
    double ans[n], *ap;
    double tmp;
    storage_type *ellp;
    int i, j, k, m;
    if (relax <= 0) {
      for (i = 0; i < n; i++)
	phi[i] = (vals[i] != 0) ? exp(t*vals[i]) : 1;
    } else {
      for (i = 0; i < n; i++)
	phi[i] = (vals[i] != 0) ? pow(1-relax*vals[i], -t/relax) : 1;
    }
    for (j = 0, ellp = ell; j < nlocus; j++, ellp += n) {
      for (i = 0, ap = ans; i < n; i++, ap++) {
	for (m = 0, *ap = 0; m < n; m++) {
	  for (k = 0, tmp = 0; k < n; k++){
	    tmp += rv[i+n*k]*phi[k]*lv[k+n*m];
	  }
	  tmp = (tmp > 0) ? tmp : 0;
	  *ap += tmp*ellp[m];
	}
      }
      for(i = 0; i < n; i++) ellp[i] = (storage_type) ans[i];
    }
  }
}


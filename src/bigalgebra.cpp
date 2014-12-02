#include <string>
#include <iostream>
#include <cmath>
#include "bigmemory/BigMatrix.h"

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#ifdef REFBLAS
#include "refblas64longlong.h"
#define INT long long
#else
#include <R_ext/BLAS.h>
#define INT int
#endif

#define BLOCKSIZE (8)

// Declare LAPACK functions
extern "C"
{
  void dgeqrf_ (int *M, int *N, double *Y, int *LDA, double *TAU,
                    double *WORK, int *LWORK, int *INFO); 
}

#ifdef __cplusplus
extern "C"
{
#endif

  double *make_double_ptr (SEXP matrix, SEXP isBigMatrix);

  SEXP dgemm_wrapper (SEXP TRANSA, SEXP TRANSB, SEXP M, SEXP N, SEXP K,
                      SEXP ALPHA, SEXP A, SEXP LDA, SEXP B, SEXP LDB,
                      SEXP BETA, SEXP C, SEXP LDC, SEXP A_isBM, SEXP B_isBM,
                      SEXP C_isBM, SEXP C_offset);
  SEXP daxpy_wrapper (SEXP N, SEXP A, SEXP X, SEXP Y, SEXP X_isBM);
  SEXP dgeqrf_wrapper (SEXP M, SEXP N, SEXP Y, SEXP LDA, SEXP TAU, SEXP WORK,
                        SEXP LWORK, SEXP INFO, SEXP A_isBM, SEXP TAU_isBM, 
                        SEXP WORK_isBM);
  SEXP dgeemm_wrapper (SEXP N, SEXP X, SEXP Y, SEXP Z, SEXP X_isBM, SEXP Y_isBM);
  SEXP dgeemd_wrapper (SEXP N, SEXP X, SEXP Y, SEXP Z, SEXP X_isBM, SEXP Y_isBM);
  SEXP dadd(SEXP N, SEXP ALPHA, SEXP Y, SEXP Y_isBM, SEXP SIGN, SEXP ALPHA_LHS);
  SEXP dgesmd_wrapper (SEXP N, SEXP A, SEXP Y, SEXP Y_isBM, SEXP ALPHA_LHS);
  
  // Generic Math functions
  SEXP dgepow(SEXP N, SEXP EXP, SEXP Y, SEXP Y_isBM);
  
  // Logarithm functions
  SEXP dgeclog(SEXP N, SEXP Y, SEXP Y_isBM);
  SEXP dgelog(SEXP N, SEXP BASE, SEXP Y, SEXP Y_isBM);
  SEXP dgeexp(SEXP N, SEXP Y, SEXP Y_isBM);
  
  // Define trignometric functions
  SEXP dgecosh(SEXP N, SEXP Y, SEXP Y_isBM);
  SEXP dgetanh(SEXP N, SEXP Y, SEXP Y_isBM);
  SEXP dgesinh(SEXP N, SEXP Y, SEXP Y_isBM);

#ifdef __cplusplus
}
#endif


/* Pointer utility, returns a double pointer for either a BigMatrix or a
 * standard R matrix.
 */
double *
make_double_ptr (SEXP matrix, SEXP isBigMatrix)
{
  double *matrix_ptr;

  if (LOGICAL_VALUE (isBigMatrix) == (Rboolean) TRUE)   // Big Matrix
    {
      SEXP address = GET_SLOT (matrix, install ("address"));
      BigMatrix *pbm =
        reinterpret_cast < BigMatrix * >(R_ExternalPtrAddr (address));
      if (!pbm)
        return (NULL);

      // Check that have acceptable big.matrix
      if (pbm->row_offset () > 0 && pbm->ncol () > 1)
        {
          std::string errMsg =
            string ("sub.big.matrix objects cannoth have row ") +
            string
            ("offset greater than zero and number of columns greater than 1");
          Rf_error (errMsg.c_str ());
          return (NULL);
        }

      index_type offset = pbm->nrow () * pbm->col_offset ();
      matrix_ptr = reinterpret_cast < double *>(pbm->matrix ()) + offset;
    }
  else                          // Regular R Matrix
    {
      matrix_ptr = NUMERIC_DATA (matrix);
    }

  return (matrix_ptr);
};


/* Wrappers for miscellaneous BLAS and LAPACK routines. */
SEXP
dgemm_wrapper (SEXP TRANSA, SEXP TRANSB, SEXP M, SEXP N, SEXP K,
               SEXP ALPHA, SEXP A, SEXP LDA, SEXP B, SEXP LDB, SEXP BETA,
               SEXP C, SEXP LDC, SEXP A_isBM, SEXP B_isBM, SEXP C_isBM,
               SEXP C_offset)
{
  long j = *(DOUBLE_DATA (C_offset));
  double *pA = make_double_ptr (A, A_isBM);
  double *pB = make_double_ptr (B, B_isBM);
  double *pC;
  SEXP ans;
  INT MM = (INT) * (DOUBLE_DATA (M));
  INT NN = (INT) * (DOUBLE_DATA (N));
  INT KK = (INT) * (DOUBLE_DATA (K));
  INT LDAA = (INT) * (DOUBLE_DATA (LDA));
  INT LDBB = (INT) * (DOUBLE_DATA (LDB));
  INT LDCC = (INT) * (DOUBLE_DATA (LDC));
  if(LOGICAL_VALUE(C_isBM) == (Rboolean) TRUE)
  {
/* Return results in a big matrix */
    pC = make_double_ptr (C, C_isBM) + j;
    PROTECT(ans = C);
  } else {
/* Allocate an output R matrix and return results there
   XXX Add check for size of MM and NN XXX 
 */
    PROTECT(ans = allocMatrix(REALSXP, (int)MM, (int)NN));
    pC = NUMERIC_DATA(ans);
  }
/* An example of an alternate C-blas interface (e.g., ACML) */
#ifdef CBLAS
  dgemm (*((char *) CHARACTER_VALUE (TRANSA)),
         *((char *) CHARACTER_VALUE (TRANSB)),
         MM, NN, KK, *(NUMERIC_DATA (ALPHA)), pA, LDAA, pB,
         LDBB, *(NUMERIC_DATA (BETA)), pC, LDCC);
#elif REFBLAS
/* Standard Fortran interface without underscoring */
  int8_dgemm ((char *) CHARACTER_VALUE (TRANSA),
         (char *) CHARACTER_VALUE (TRANSB),
         &MM, &NN, &KK, NUMERIC_DATA (ALPHA), pA, &LDAA, pB,
         &LDBB, NUMERIC_DATA (BETA), pC, &LDCC);
#else
/* Standard Fortran interface from R's blas */
  dgemm_ ((char *) CHARACTER_VALUE (TRANSA),
         (char *) CHARACTER_VALUE (TRANSB),
         &MM, &NN, &KK, NUMERIC_DATA (ALPHA), pA, &LDAA, pB,
         &LDBB, NUMERIC_DATA (BETA), pC, &LDCC);
#endif
  unprotect(1);
  return ans;
}



/* Compute A*X + Y for scalar a, vectors X and Y of length N.
 * Y must be a big.matrix, X can be an R vector or big.matrix.
 * The contents of Y are *replaced* by this routine and a reference
 * to Y is returned.
 */
SEXP
daxpy_wrapper (SEXP N, SEXP A, SEXP X, SEXP Y, SEXP X_isBM)
{
  SEXP ans, Tr;
  double *pY;
  double *pA = DOUBLE_DATA(A);
  double *pX = make_double_ptr (X, X_isBM);
  INT incx = 1;
  INT incy = 1;
  INT NN = (INT) * (DOUBLE_DATA (N));
  PROTECT(ans = Y);
  PROTECT(Tr = allocVector(LGLSXP, 1));
  LOGICAL(Tr)[0] = 1;
  pY = make_double_ptr (Y, Tr);
/* An example of an alternate C-blas interface (e.g., ACML) */
#ifdef CBLAS
  daxpy_ (NN, pA, pX, incx, pY, incy);
#elif REFBLAS
/* Standard Fortran interface without underscoring */
  int8_daxpy (&NN, pA, pX, &incx, pY, &incy);
#else
/* Standard Fortran interface from R's blas */
  daxpy_ (&NN, pA, pX, &incx, pY, &incy);
#endif
  unprotect(2);
  return ans;
}

// Note that the following two functions should be updated to use
// the proper BLAS functions.

/* Note, I am not certain that this function is necessary.
 * Currently, modified to have a scalar converted to a big.matrix
 * in order to use daxpy
 */


SEXP
dadd(SEXP N, SEXP ALPHA, SEXP Y, SEXP Y_isBM, SEXP SIGN, SEXP ALPHA_LHS) {
  SEXP ans;
  double *pY;
  INT NN = (INT) * (DOUBLE_DATA(N));
  double alpha = *(DOUBLE_DATA(ALPHA));
  int alpha_lhs = *(INTEGER_DATA(ALPHA_LHS));
  double sign = *(DOUBLE_DATA(SIGN));
  pY = make_double_ptr(Y, Y_isBM);
  PROTECT(ans = Y);
  if (alpha_lhs) 
  {
    for (INT i=0; i < NN; ++i)
    {
      pY[i] = alpha + sign * pY[i];
    }
  }
  else
  {
    for (INT i=0; i < NN; ++i)
    {
      pY[i] += sign*alpha;
    }
  }
  unprotect(1);
  return ans;
}




/* Compute the QR decomposition of a big.matrix
 * Y must be a big.matrix, X can be an R vector or big.matrix.
 * The contents of Y are *replaced* by this routine and a reference
 * to Y is returned.
 */
SEXP
dgeqrf_wrapper (SEXP M, SEXP N, SEXP Y, SEXP LDA, SEXP TAU, SEXP WORK,
                SEXP LWORK, SEXP INFO, SEXP A_isBM, SEXP TAU_isBM, 
                SEXP WORK_isBM)
{
  SEXP ans, Tr;
  double *pY;
  double *pTAU = make_double_ptr (TAU, TAU_isBM);
  double *pWORK = make_double_ptr (WORK, WORK_isBM);
  INT NN = (INT) * (DOUBLE_DATA (N));
  INT MM = (INT) * (DOUBLE_DATA (M));
  INT LDAi = (INT) * (DOUBLE_DATA (LDA));
  INT LWORKi = (INT) * (DOUBLE_DATA (LWORK));
  INT INFOi = (INT) * (DOUBLE_DATA (INFO));
  
  PROTECT(ans = Y);
  PROTECT(Tr = allocVector(LGLSXP, 1));
  LOGICAL(Tr)[0] = 1;
  pY = make_double_ptr (Y, Tr);
  
/* C-blas and REFBLAS untested */
/* An example of an alternate C-blas interface (e.g., ACML) */
#ifdef CBLAS
  dgeqrf (MM, NN, pY, LDAi, pTAU, pWORK, LWORKi, INFOi);
#elif REFBLAS
/* Standard Fortran interface without underscoring */
  int8_dgeqrf (&MM, &NN, pY, &LDAi, pTAU, pWORK, &LWORKi, &INFOi);
#else
/* Connect to LAPACK Fortran interface
 * The function is declared at start of this file
 */
  dgeqrf_ (&MM, &NN, pY, &LDAi, pTAU, pWORK, &LWORKi, &INFOi);
#endif
  unprotect(2);
  return ans;
}


// Unsure if there is a 'proper' BLAS function for the following functions

// element-wise matrix multiplcation
SEXP
dgeemm_wrapper (SEXP N, SEXP X, SEXP Y, SEXP Z, SEXP X_isBM, SEXP Y_isBM){
  SEXP ans, Tr;
  double *pZ;
  double *pY = make_double_ptr (Y, Y_isBM);
  double *pX = make_double_ptr (X, X_isBM);
  unsigned int NN = (unsigned int) * (DOUBLE_DATA(N));
//  pY = make_double_ptr(Y, Y_isBM);
  
  PROTECT(ans = Z);
  PROTECT(Tr = allocVector(LGLSXP, 1));
  LOGICAL(Tr)[0] = 1;
  pZ = make_double_ptr (Z, Tr);
  
  unsigned int i = 0;
  unsigned int blocklimit;
  
  /* The limit may not be divisible by BLOCKSIZE, 
   * go as near as we can first, then tidy up.
   */ 
  blocklimit = ( NN / BLOCKSIZE ) * BLOCKSIZE;
  
  while( i < blocklimit )
  {
    pZ[i] = pX[i]*pY[i];
    pZ[i+1] = pX[i+1]*pY[i+1];
    pZ[i+2] = pX[i+2]*pY[i+2];
    pZ[i+3] = pX[i+3]*pY[i+3];
    pZ[i+4] = pX[i+4]*pY[i+4];
    pZ[i+5] = pX[i+5]*pY[i+5];
    pZ[i+6] = pX[i+6]*pY[i+6];
    pZ[i+7] = pX[i+7]*pY[i+7];
    
    // update counter
    i+=8;
  }
  
  // finish remaining elements
  if( i < NN ) 
    { 
        /* Jump into the case at the place that will allow
         * us to finish off the appropriate number of items. 
         */ 

        switch( NN - i ) 
        { 
            case 7 : pZ[i] = pX[i]*pY[i]; i++; 
            case 6 : pZ[i] = pX[i]*pY[i]; i++; 
            case 5 : pZ[i] = pX[i]*pY[i]; i++; 
            case 4 : pZ[i] = pX[i]*pY[i]; i++; 
            case 3 : pZ[i] = pX[i]*pY[i]; i++; 
            case 2 : pZ[i] = pX[i]*pY[i]; i++; 
            case 1 : pZ[i] = pX[i]*pY[i]; 
        }
    } 
    unprotect(2);
    return ans;
}


// element-wise matrix division
SEXP
dgeemd_wrapper (SEXP N, SEXP X, SEXP Y, SEXP Z, SEXP X_isBM, SEXP Y_isBM){
  SEXP ans, Tr;
  double *pZ;
  double *pY = make_double_ptr (Y, Y_isBM);
  double *pX = make_double_ptr (X, X_isBM);
  unsigned int NN = (unsigned int) * (DOUBLE_DATA(N));
//  pY = make_double_ptr(Y, Y_isBM);
  
  PROTECT(ans = Z);
  PROTECT(Tr = allocVector(LGLSXP, 1));
  LOGICAL(Tr)[0] = 1;
  pZ = make_double_ptr (Z, Tr);
  
  unsigned int i = 0;
  //unsigned int limit = 33;
  unsigned int blocklimit;
  
  /* The limit may not be divisible by BLOCKSIZE, 
   * go as near as we can first, then tidy up.
   */ 
  blocklimit = ( NN / BLOCKSIZE ) * BLOCKSIZE;
  
  while( i < blocklimit )
  {
    pZ[i] = pX[i] / pY[i];
    pZ[i+1] = pX[i+1] / pY[i+1];
    pZ[i+2] = pX[i+2] / pY[i+2];
    pZ[i+3] = pX[i+3] / pY[i+3];
    pZ[i+4] = pX[i+4] / pY[i+4];
    pZ[i+5] = pX[i+5] / pY[i+5];
    pZ[i+6] = pX[i+6] / pY[i+6];
    pZ[i+7] = pX[i+7] / pY[i+7];
    
    // update counter
    i+=8;
  }
  
  // finish remaining elements
  if( i < NN ) 
    { 
        /* Jump into the case at the place that will allow
         * us to finish off the appropriate number of items. 
         */ 

        switch( NN - i ) 
        { 
            case 7 : pZ[i] = pX[i] / pY[i]; i++; 
            case 6 : pZ[i] = pX[i] / pY[i]; i++; 
            case 5 : pZ[i] = pX[i] / pY[i]; i++; 
            case 4 : pZ[i] = pX[i] / pY[i]; i++; 
            case 3 : pZ[i] = pX[i] / pY[i]; i++; 
            case 2 : pZ[i] = pX[i] / pY[i]; i++; 
            case 1 : pZ[i] = pX[i] / pY[i]; 
        }
    } 
    unprotect(2);
    return ans;
}


// Scalar-matrix division
SEXP
dgesmd_wrapper (SEXP N, SEXP A, SEXP Y, SEXP Y_isBM, SEXP ALPHA_LHS) {
  SEXP ans;
  double *pY;
  INT NN = (INT) * (DOUBLE_DATA(N));
  INT ALPHA = (INT) * (DOUBLE_DATA(A));
  int alpha_lhs = *(INTEGER_DATA(ALPHA_LHS));
  pY = make_double_ptr(Y, Y_isBM);
  PROTECT(ans = Y);
 
  unsigned int i = 0;
  unsigned int blocklimit;
  blocklimit = ( NN / BLOCKSIZE ) * BLOCKSIZE;
  
  if(alpha_lhs == 0){
    while( i < blocklimit )
    {
      pY[i] = pY[i]/ALPHA;
      pY[i+1] = pY[i+1]/ALPHA;
      pY[i+2] = pY[i+2]/ALPHA;
      pY[i+3] = pY[i+3]/ALPHA;
      pY[i+4] = pY[i+4]/ALPHA;
      pY[i+5] = pY[i+5]/ALPHA;
      pY[i+6] = pY[i+6]/ALPHA;
      pY[i+7] = pY[i+7]/ALPHA;
      
      // update counter
      i+=8;
    }
    
    // finish remaining elements
    if( i < NN ) 
      { 
          switch( NN - i ) 
          { 
              case 7 : pY[i] = pY[i] / ALPHA; i++; 
              case 6 : pY[i] = pY[i] / ALPHA; i++; 
              case 5 : pY[i] = pY[i] / ALPHA; i++; 
              case 4 : pY[i] = pY[i] / ALPHA; i++; 
              case 3 : pY[i] = pY[i] / ALPHA; i++; 
              case 2 : pY[i] = pY[i] / ALPHA; i++; 
              case 1 : pY[i] = pY[i] / ALPHA; 
          }
      } 
  }
  else
  {
    while( i < blocklimit )
    {
      pY[i] = ALPHA / pY[i];
      pY[i+1] = ALPHA / pY[i+1];
      pY[i+2] = ALPHA / pY[i+2];
      pY[i+3] = ALPHA / pY[i+3];
      pY[i+4] = ALPHA / pY[i+4];
      pY[i+5] = ALPHA / pY[i+5];
      pY[i+6] = ALPHA / pY[i+6];
      pY[i+7] = ALPHA / pY[i+7];
      
      // update counter
      i+=8;
    }
    
    // finish remaining elements
    if( i < NN ) 
      { 
          switch( NN - i ) 
          { 
              case 7 : pY[i] = ALPHA / pY[i]; i++; 
              case 6 : pY[i] = ALPHA / pY[i]; i++; 
              case 5 : pY[i] = ALPHA / pY[i]; i++; 
              case 4 : pY[i] = ALPHA / pY[i]; i++; 
              case 3 : pY[i] = ALPHA / pY[i]; i++; 
              case 2 : pY[i] = ALPHA / pY[i]; i++; 
              case 1 : pY[i] = ALPHA / pY[i]; 
          }
      }
  }
  
  unprotect(1);
  return ans;
}


// Power
SEXP
dgepow(SEXP N, SEXP EXP, SEXP Y, SEXP Y_isBM) {
  SEXP ans;
  double *pY;
  INT NN = (INT) * (DOUBLE_DATA(N));
  double PP = (double) * (DOUBLE_DATA (EXP));
  pY = make_double_ptr(Y, Y_isBM);
  PROTECT(ans = Y);
 
  for (INT i=0; i < NN; ++i)
  {
    pY[i] = pow(pY[i], PP);
  }
  
  unprotect(1);
  return ans;
}


// common logarithm
SEXP
dgeclog(SEXP N, SEXP Y, SEXP Y_isBM) {
  SEXP ans;
  double *pY;
  INT NN = (INT) * (DOUBLE_DATA(N));
  pY = make_double_ptr(Y, Y_isBM);
  PROTECT(ans = Y);
 
  for (INT i=0; i < NN; ++i)
  {
    pY[i] = log10(pY[i]);
  }
  
  unprotect(1);
  return ans;
}

// base logarithm
SEXP
dgelog(SEXP N, SEXP BASE, SEXP Y, SEXP Y_isBM) {
  SEXP ans;
  double *pY;
  INT NN = (INT) * (DOUBLE_DATA(N));
  double BB = (double) * (DOUBLE_DATA (BASE));
  pY = make_double_ptr(Y, Y_isBM);
  PROTECT(ans = Y);
 
  for (INT i=0; i < NN; ++i)
  {
    pY[i] = log10(pY[i])/log10(BB);
  }
  
  unprotect(1);
  return ans;
}

// Exponential function
SEXP
dgeexp(SEXP N, SEXP Y, SEXP Y_isBM) {
  SEXP ans;
  double *pY;
  INT NN = (INT) * (DOUBLE_DATA(N));
  pY = make_double_ptr(Y, Y_isBM);
  PROTECT(ans = Y);
 
  for (INT i=0; i < NN; ++i)
  {
    pY[i] = exp(pY[i]);
  }
  
  unprotect(1);
  return ans;
}

// hyperbolic tangent
SEXP
dgetanh(SEXP N, SEXP Y, SEXP Y_isBM) {
  SEXP ans;
  double *pY;
  INT NN = (INT) * (DOUBLE_DATA(N));
  pY = make_double_ptr(Y, Y_isBM);
  PROTECT(ans = Y);
 
  for (INT i=0; i < NN; ++i)
  {
    pY[i] = tanh(pY[i]);
  }
  
  unprotect(1);
  return ans;
}

// hyperbolic cosine
SEXP
dgecosh(SEXP N, SEXP Y, SEXP Y_isBM) {
  SEXP ans;
  double *pY;
  INT NN = (INT) * (DOUBLE_DATA(N));
  pY = make_double_ptr(Y, Y_isBM);
  PROTECT(ans = Y);
 
  for (INT i=0; i < NN; ++i)
  {
    pY[i] = cosh(pY[i]);
  }
  
  unprotect(1);
  return ans;
}

// hyperbolic sine
SEXP
dgesinh(SEXP N, SEXP Y, SEXP Y_isBM) {
  SEXP ans;
  double *pY;
  INT NN = (INT) * (DOUBLE_DATA(N));
  pY = make_double_ptr(Y, Y_isBM);
  PROTECT(ans = Y);
 
  for (INT i=0; i < NN; ++i)
  {
    pY[i] = sinh(pY[i]);
  }
  
  unprotect(1);
  return ans;
}
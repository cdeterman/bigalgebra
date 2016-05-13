// prevent R include files from defining length
#define R_NO_REMAP

// add BigMatrix 
#include "bigmemory/BigMatrix.h"
#include "bigmemory/MatrixAccessor.hpp"

// condition for 64-bit reference BLAS
#ifdef REFBLAS
#include "refblas64longlong.h"
#define INT long long
#else
// no longer required with RcppArmadillo
//#include <R_ext/BLAS.h>
//#include <R_ext/Lapack.h>
#define INT int
#endif


#define ARMA_NO_DEBUG

/* additional macros that could potentially
 * help solve the 64-bit issue that requires
 * the REFBLAS component
 */
 
//#define ARMA_64BIT_WORD
//#define ARMA_BLAS_LONG_LONG

#include <Rdefines.h>

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

#include "bigmatrix_arma.h"

using namespace Rcpp;

#ifdef __cplusplus
extern "C"
{
#endif

  double *make_double_ptr (SEXP matrix, SEXP isBigMatrix);

#ifdef __cplusplus
}
#endif

#if defined(CBLAS) || defined(REFBLAS)
#define NON_R_BLAS=1
#endif

/* The required roxygen elements for the NAMESPACE to
 * be created correctly.
 */
 
//' @useDynLib bigalgebra
//' @importFrom Rcpp evalCpp


/* Pointer utility, returns a double pointer for either a BigMatrix or a
 * standard R matrix.
 */
double *
make_double_ptr (SEXP matrix, SEXP isBigMatrix)
{
  double *matrix_ptr;

  if (Rf_asLogical (isBigMatrix) == (Rboolean) TRUE)   // Big Matrix
    {
      SEXP address = R_do_slot (matrix, Rf_install ("address"));
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
}


/* Wrappers for miscellaneous BLAS and LAPACK routines. */

/* I have restructured the following functions as the RcppArmadillo
 * implementation requires far less variables to apply.
 * As such, the conditional compilation has been altered to
 * declare the appropriate variables depending upon the
 * current BLAS implementation used
 */
 
// [[Rcpp::export]]
SEXP
dgemm_wrapper (SEXP TRANSA, SEXP TRANSB, SEXP M, SEXP N, SEXP K,
               SEXP ALPHA, SEXP A, SEXP LDA, SEXP B, SEXP LDB, SEXP BETA,
               SEXP C, SEXP LDC, bool A_isBM, bool B_isBM, bool C_isBM,
               SEXP C_offset)
{
  
  #ifdef NON_R_BLAS
    SEXP ans;
    long j = *(DOUBLE_DATA (C_offset));
    double *pA = make_double_ptr (A, A_isBM);
    double *pB = make_double_ptr (B, B_isBM); 
    double *pC;
    INT KK = (INT) * (DOUBLE_DATA (K));
    INT LDAA = (INT) * (DOUBLE_DATA (LDA));
    INT LDBB = (INT) * (DOUBLE_DATA (LDB));
    INT LDCC = (INT) * (DOUBLE_DATA (LDC));
    INT MM = (INT) * (DOUBLE_DATA (M));
    INT NN = (INT) * (DOUBLE_DATA (N));
  
    if(Rf_asLogical(C_isBM) == (Rboolean) TRUE)
      {
        /* Return results in a big matrix */
        pC = make_double_ptr (C, C_isBM) + j;
        PROTECT(ans = C);
      } else {
        /* Allocate an output R matrix and return results there
         * XXX Add check for size of MM and NN XXX 
         */
        PROTECT(ans = Rf_allocMatrix(REALSXP, (int)MM, (int)NN));
        pC = NUMERIC_DATA(ans);
      }
      
    /* An example of an alternate C-blas interface (e.g., ACML) */
    #ifdef CBLAS
    dgemm (*((char *) CHAR(Rf_asChar (TRANSA))),
           *((char *) CHAR(Rf_asChar (TRANSB))),
           MM, NN, KK, *(NUMERIC_DATA (ALPHA)), pA, LDAA, pB,
           LDBB, *(NUMERIC_DATA (BETA)), pC, LDCC);
    #elif REFBLAS
      /* Standard Fortran interface without underscoring */
      int8_dgemm ((char *) CHAR(Rf_asChar (TRANSA)),
             (char *) CHAR(Rf_asChar (TRANSB)),
             &MM, &NN, &KK, NUMERIC_DATA (ALPHA), pA, &LDAA, pB,
             &LDBB, NUMERIC_DATA (BETA), pC, &LDCC);
    #else
    #error "BLAS format not supported"
    #endif
    
    Rf_unprotect(1);
    return ans;
    
  #else

    const arma::mat Am( A_isBM ? ConvertBMtoArma(A) : as<arma::mat>(A) );
    const arma::mat Bm( B_isBM ? ConvertBMtoArma(B) : as<arma::mat>(B) );
    arma::mat Cm( C_isBM ? ConvertBMtoArma(C) : as<arma::mat>(C) );
    
    /* RcppArmadillo BLAS interface */
    // RcppArmadillo Matrix Multiplication
    Cm = Am * Bm;
        
    if(C_isBM)
    {
      return C;
    }else{
      // return normal R matrix
      return Rcpp::wrap(Cm);
    } 
  #endif
}


/* Compute A*X + Y for scalar a, vectors X and Y of length N.
 * Y must be a big.matrix, X can be an R vector or big.matrix.
 * The contents of Y are *replaced* by this routine and a reference
 * to Y is returned.
 */
 // [[Rcpp::export]]
SEXP
daxpy_wrapper (SEXP N, SEXP A, SEXP X, SEXP Y, bool X_isBM)
{
  
  #ifdef NON_R_BLAS
    SEXP ans, Tr;
    double *pY;
    double *pA = DOUBLE_DATA(A);
    double *pX = make_double_ptr (X, X_isBM);
    INT incx = 1;
    INT incy = 1;
    INT NN = (INT) * (DOUBLE_DATA (N));
    PROTECT(ans = Y);
    PROTECT(Tr = Rf_allocVector(LGLSXP, 1));
    LOGICAL(Tr)[0] = 1;
    pY = make_double_ptr (Y, Tr);
    /* An example of an alternate C-blas interface (e.g., ACML) */
    #ifdef CBLAS
      daxpy_ (NN, pA, pX, incx, pY, incy);
    #elif REFBLAS
    /* Standard Fortran interface without underscoring */
      int8_daxpy (&NN, pA, pX, &incx, pY, &incy);
    #else
      #error "BLAS format not supported"
    #endif
    
    Rf_unprotect(2);
    return ans;
    
  #else 
    /* RcppArmadillo 'BLAS' implementation */
    
    // convert to arma matrices
    const arma::mat Xm( X_isBM ? ConvertBMtoArma(X) : as<arma::mat>(X) );
    arma::mat Ym( X_isBM ? ConvertBMtoArma(Y) : as<arma::mat>(Y) );
    
    double ALPHA = Rcpp::as<double>(A);
    
    // RcppArmadillo Matrix Addition or Subtraction
    Ym = ALPHA * Xm + Ym;
    
    return Y;
   
  #endif
}
  
// [[Rcpp::export]]
SEXP dpotrf_wrapper(SEXP UPLO, SEXP N, SEXP A, SEXP LDA, SEXP INFO, bool A_isBM)
{
#ifdef NON_R_BLAS
  SEXP ans;
  const char *_UPLO = CHAR(Rf_asChar(UPLO));
  INT _N = (INT)* (DOUBLE_DATA(N));
  double *_A = make_double_ptr(A, A_isBM);
  INT _LDA = (INT) *(DOUBLE_DATA(LDA));
  INT _INFO = (INT) *(DOUBLE_DATA(INFO));
  /* An example of an alternate C-blas interface (e.g., ACML) */
  #ifdef CBLAS
    dpotrf (_UPLO, &_N, _A, &_LDA, &_INFO);
  #elif REFBLAS
    /* Standard Fortran interface without underscoring */
    int8_dpotrf (_UPLO, &_N, _A, &_LDA, &_INFO);
  #else
    #error "BLAS format not supported"
  #endif
    PROTECT(ans = A);
    Rf_unprotect(1);
    return ans;
  #else
    
    // convert to armadillo matrix
    arma::mat Am( A_isBM ? ConvertBMtoArma(A) : as<arma::mat>(A) );
    
    Am = arma::chol(Am);
    return A;
  #endif
}

// Note that the following two functions should be updated to use
// the proper BLAS functions.

/* Note, I am not certain that this function is necessary.
 * Currently, modified to have a scalar converted to a big.matrix
 * in order to use daxpy
 */
 
/* This has been updated to utilize RcppArmadillo which
 * provides an interface to BLAS and LAPACK functions
 */

// [[Rcpp::export]]
SEXP
dadd_wrapper(SEXP ALPHA, SEXP Y, bool Y_isBM, SEXP SIGN) {
    // convert to arma matrices
    arma::mat Ym( Y_isBM ? ConvertBMtoArma(Y) : as<arma::mat>(Y) );


    double SCALAR = Rcpp::as<double>(ALPHA);
    int sign = Rcpp::as<int>(SIGN);
    
    // RcppArmadillo Matrix Addition or Subtraction
    Ym = sign * SCALAR + Ym;
    
    return Y;
}


/* Compute the QR decomposition of a big.matrix
 * Y, Q, & R must be of class big.matrix
 * 
 * Also, the output is currently different between armadillo and
 * direct LAPACK.  Armadillo output is consistent with other
 * functions such as numpy.linalg.qr or scipy.linalg.qr in python.
 * As such, I currently am leaning towards this output as opposed
 * to the 'compact' form returned by the base R 'qr' function.
 */

// [[Rcpp::export]]
SEXP
dgeqrf_wrapper (SEXP Y, SEXP Q, SEXP R)
{
    const arma::mat Ym = ConvertBMtoArma(Y);
    arma::mat Qm = ConvertBMtoArma(Q);
    arma::mat Rm = ConvertBMtoArma(R);
    
    // QR decomposition
    /* arma::qr returns a bool set to false if QR fails
     * need to have an exception catch.
     * NOTE - This is untested!!! Need a matrix that will fail
     */
    if(!arma::qr(Qm,Rm,Ym)){
      Rcpp::stop("QR decomposition failed");
    };
    
    // return list of Q and R
    return Rcpp::List::create(Named("Q") = Q, Named("R") = R);
}


/* No specific BLAS functions for the following so all
 * are strictly armadillo.
 */

// element-wise matrix multiplcation
// [[Rcpp::export]]
SEXP
dgeemm_wrapper (SEXP X, SEXP Y, bool X_isBM){
  
    // convert to arma matrices
    const arma::mat Xm( X_isBM ? ConvertBMtoArma(X) : as<arma::mat>(X) );

    // Y will always be a big.matrix
    arma::mat Ym = ConvertBMtoArma(Y);
    
    // RcppArmadillo Element-wise Matrix Multiplication
    Ym = Xm % Ym;
    
    return Y;
}


// element-wise matrix division
// [[Rcpp::export]]
SEXP
dgeemd_wrapper (SEXP X, SEXP Y, bool X_isBM){
    // convert to arma matrices
    const arma::mat Xm( X_isBM ? ConvertBMtoArma(X) : as<arma::mat>(X) );

    // Y will always be a big.matrix b/c created by anon_matrix
//    arma::mat Ym( Y_isBM ? ConvertBMtoArma(Y) : as<arma::mat>(Y) );
    arma::mat Ym = ConvertBMtoArma(Y);
    
    // RcppArmadillo Element-wise Matrix Division
    Ym = Xm / Ym;
    
    return Y;
}


// Scalar-matrix division
// [[Rcpp::export]]
SEXP
dgesmd_wrapper (SEXP A, SEXP Y, bool Y_isBM, int ALPHA_LHS) {

    INT ALPHA = (INT) * (DOUBLE_DATA(A));
  
    // Y will always be a big.matrix
    arma::mat Ym = ConvertBMtoArma(Y);
    
    // RcppArmadillo Element-wise Scalar Division
    if(ALPHA_LHS){
        Ym = ALPHA / Ym;
    }else{      
        Ym = Ym / ALPHA;
    }
    
    return Y;
}


//' @name t_wrapper
//' @title Big Matrix Transposition
//' @param X big.matrix object
//' @param Y big.matrix object (to be filled)
//' @param LOW_MEM boolean option to choose slower, low memory option
// [[Rcpp::export]]
SEXP
t_wrapper (SEXP X, SEXP Y){
    const arma::mat Xm = ConvertBMtoArma(X);
    arma::mat Ym = ConvertBMtoArma(Y);
    
    // matrix transposition
    Ym = trans(Xm);

    return Y;
}

//' @name t_inplace_wrapper
//' @title Big Matrix In Place Transposition
//' @param X big.matrix object
//' @param LOW_MEM boolean option to choose slower, low memory option
// [[Rcpp::export]]
SEXP
t_inplace_wrapper (SEXP X, SEXP LOW_MEM){

    arma::mat Xm = ConvertBMtoArma(X);
    
    if(Rf_asLogical(LOW_MEM) == (Rboolean) TRUE){
        arma::inplace_trans(Xm);
    }else{
        arma::inplace_trans(Xm, "lowmem");
    }
    
    return X;
}

// Power
// [[Rcpp::export]]
SEXP
dgepow_wrapper(SEXP EXP, SEXP Y) {
    // Y will always be a big.matrix
    arma::mat Ym = ConvertBMtoArma(Y);
      
    // RcppArmadillo Element-wise Scalar Division
    Ym = pow(Ym, Rcpp::as<double>(EXP));
    
    return Y;
}


// common logarithm
// [[Rcpp::export]]
SEXP
dgeclog_wrapper(SEXP Y) {
    // Y will always be a big.matrix
    arma::mat Ym = ConvertBMtoArma(Y);
    
    // RcppArmadillo Element-wise Scalar Division
    Ym = log10(Ym);
    
    return Y;
}

// base logarithm
// [[Rcpp::export]]
SEXP
dgelog_wrapper(SEXP BASE, SEXP Y) {
    // Y will always be a big.matrix
    arma::mat Ym = ConvertBMtoArma(Y);
    
    // RcppArmadillo Element-wise Scalar Division
    Ym = log10(Ym)/log10(Rcpp::as<double>(BASE));
    
    return Y;
}

// Exponential function
// [[Rcpp::export]]
SEXP
dgeexp_wrapper(SEXP Y) {
    // Y will always be a big.matrix
    arma::mat Ym = ConvertBMtoArma(Y);
    
    // RcppArmadillo Element-wise Scalar Division
    Ym = exp(Ym);
    
    return Y;
}

// hyperbolic tangent
// [[Rcpp::export]]
SEXP
dgetanh_wrapper(SEXP Y) {
    // Y will always be a big.matrix
    arma::mat Ym = ConvertBMtoArma(Y);
    
    // RcppArmadillo Element-wise Scalar Division
    Ym = tanh(Ym);
    
    return Y;
}

// hyperbolic cosine
// [[Rcpp::export]]
SEXP
dgecosh_wrapper(SEXP Y) {
    // Y will always be a big.matrix
    arma::mat Ym = ConvertBMtoArma(Y);
    
    // RcppArmadillo Element-wise Scalar Division
    Ym = cosh(Ym);
  
    return Y;
}

// hyperbolic sine
// [[Rcpp::export]]
SEXP
dgesinh_wrapper(SEXP Y) {
    // Y will always be a big.matrix
    arma::mat Ym = ConvertBMtoArma(Y);
    
    // RcppArmadillo Element-wise Scalar Division
    Ym = sinh(Ym);
    
    return Y;
}

// eigen values
// [[Rcpp::export]]
Rcpp::List
eigen_wrapper(SEXP X, SEXP EIG_VECS, SEXP only_values){
  
  arma::mat Xm = ConvertBMtoArma(X);
  
  if(Rf_asLogical(only_values) == (Rboolean) TRUE){
      
      arma::vec eigval;
      arma::eig_sym(eigval, Xm);
      
      // To match R return in reverse order
      // may remove this      
      eigval = arma::sort(eigval, "descend");
      
      return Rcpp::List::create(
          Named("values") = Rcpp::NumericVector(eigval.begin(), eigval.end()));
  }else{
      arma::vec eigval;
      arma::mat eigvec;
      arma::uword MM = Xm.n_cols;
      
      arma::eig_sym(eigval, eigvec, Xm);
      arma::uvec idx = arma::sort_index(eigval, "descend");
      
      // To match R return in descending order
      // may remove this
      eigval = arma::sort(eigval, "descend");
      
      arma::mat eigvec_out = ConvertBMtoArma(EIG_VECS);
      
      // this loop uses the index for the eigval sort to 
      // transfer the matrix in that sorted order
      arma::uword i = 0;
      for(arma::uvec::iterator it = idx.begin(); it != idx.end(); ++it){
        std::memcpy(eigvec_out.colptr(i), eigvec.colptr(*it), MM*sizeof(double));
        i++;
      }
      
      return Rcpp::List::create(Named("values") = Rcpp::NumericVector(eigval.begin(), eigval.end()), 
                                Named("vectors") = EIG_VECS);
  }
}

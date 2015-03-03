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
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#define INT int
#endif


#define ARMA_NO_DEBUG

/* additional macros that could potentially
 * help solve the 64-bit issue that requires
 * the REFBLAS component
 */
 
//#define ARMA_64BIT_WORD
//#define ARMA_BLAS_LONG_LONG

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

#define BLOCKSIZE (8)

using namespace Rcpp;

// Declare LAPACK functions
extern "C"
{
  //void dgeqrf_ (int *M, int *N, double *Y, int *LDA, double *TAU,
  //                  double *WORK, int *LWORK, int *INFO); 
}

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

//' @useDynLib bigalgebra
//' @importFrom Rcpp evalCpp
//' @export
// [[Rcpp::export]]
void hello(){
  std::cout << "hello" << std::endl;
}

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
 
// simple function to get access to big.matrix objects from R
Rcpp::XPtr<BigMatrix> BigMatrixXPtr(SEXP A){
    // declare as S4 object
    Rcpp::S4 As4(A);
    // pull address slot
    SEXP A_address = As4.slot("address");
    // declare as external pointer
    Rcpp::XPtr<BigMatrix> xpA(A_address);
    return xpA;
 }
       
 
arma::mat ConvertToArma(SEXP A, SEXP isBM)
{
  if(Rf_asLogical(isBM) == (Rboolean) TRUE)
      { 
        Rcpp::XPtr<BigMatrix> xpA = BigMatrixXPtr(A);
        
        arma::mat Am = arma::mat( (double*) xpA->matrix(),
                    xpA->nrow(),
                    xpA->ncol(),
                    false);
        return Am;
      }else{
        arma::mat Am = Rcpp::as<arma::mat>(A);
        return Am;
      }
}

arma::mat ConvertBMtoArma(SEXP A)
{
  Rcpp::XPtr<BigMatrix> xpA = BigMatrixXPtr(A);
        
  arma::mat Am = arma::mat( (double*) xpA->matrix(),
              xpA->nrow(),
              xpA->ncol(),
              false);
  return Am;
}
 
// [[Rcpp::export]]
SEXP
dgemm_wrapper (SEXP TRANSA, SEXP TRANSB, SEXP M, SEXP N, SEXP K,
               SEXP ALPHA, SEXP A, SEXP LDA, SEXP B, SEXP LDB, SEXP BETA,
               SEXP C, SEXP LDC, SEXP A_isBM, SEXP B_isBM, SEXP C_isBM,
               SEXP C_offset)
{

  INT MM = (INT) * (DOUBLE_DATA (M));
  INT NN = (INT) * (DOUBLE_DATA (N));
  
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
    /* RcppArmadillo BLAS interface */
    const arma::mat Am = ConvertToArma(A, A_isBM);
    const arma::mat Bm = ConvertToArma(B, B_isBM);
    
    // RcppArmadillo Matrix Multiplication
    const arma::mat Cm = Am * Bm;
    
    /* Ideally would like to avoid the following loop
     * would need to consider passing another check for
     * a separated big.matrix.  The following is required
     * for a separated big.matrix so it works in both scenarios
     */
    
    if(Rf_asLogical(C_isBM) == (Rboolean) TRUE)
    {
      // get BigMatrix external pointer
      Rcpp::XPtr<BigMatrix> xpC = BigMatrixXPtr(C);
      // create MatrixAccessor object to imput values
      MatrixAccessor<double> C_BM(*xpC);
      
      // memcpy elements to the big.matrix object
      for (int i = 0; i < MM; i++){
        std::memcpy(C_BM[i], Cm.colptr(i), NN*sizeof(double));
      }
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
daxpy_wrapper (SEXP N, SEXP A, SEXP X, SEXP Y, SEXP X_isBM)
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
    const arma::mat Xm = ConvertToArma(X, X_isBM);
    arma::mat Ym = ConvertToArma(Y, X_isBM);
    
    arma::uword MM = Xm.n_cols;
    arma::uword NN = Xm.n_rows;

    double ALPHA = Rcpp::as<double>(A);
    
    // RcppArmadillo Matrix Addition or Subtraction
    Ym = ALPHA * Xm + Ym;
    
    /* Ideally would like to avoid the following loop
     * would need to consider passing another check for
     * a separated big.matrix.  The following is required
     * for a separated big.matrix so it works in both scenarios
     */
     
    // get BigMatrix external pointer
    Rcpp::XPtr<BigMatrix> xpY = BigMatrixXPtr(Y);
    // create MatrixAccessor object to imput values
    MatrixAccessor<double> Y_BM(*xpY);
    
    // memcpy elements to the big.matrix object
    for (arma::uword i = 0; i < MM; i++){
      std::memcpy(Y_BM[i], Ym.colptr(i), NN*sizeof(double));
    }
    return Y;
   
  #endif
}
  
// [[Rcpp::export]]
SEXP dpotrf_wrapper(SEXP UPLO, SEXP N, SEXP A, SEXP LDA, SEXP INFO, SEXP A_isBM)
{
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
/* Standard Fortran interface from R's lapack */
  dpotrf_ (_UPLO, &_N, _A, &_LDA, &_INFO);
#endif
  PROTECT(ans = A);
  Rf_unprotect(1);
  return ans;
}

// Note that the following two functions should be updated to use
// the proper BLAS functions.

/* Note, I am not certain that this function is necessary.
 * Currently, modified to have a scalar converted to a big.matrix
 * in order to use daxpy
 */

// [[Rcpp::export]]
SEXP
dadd_wrapper(SEXP ALPHA, SEXP Y, SEXP Y_isBM, SEXP SIGN) {
    // convert to arma matrices
    arma::mat Ym = ConvertToArma(Y, Y_isBM);
    
    arma::uword MM = Ym.n_cols;
    arma::uword NN = Ym.n_rows;

    double SCALAR = Rcpp::as<double>(ALPHA);
    double sign = Rcpp::as<double>(SIGN);
    
    // RcppArmadillo Matrix Addition or Subtraction
    Ym = sign * SCALAR + Ym;
    
    /* Ideally would like to avoid the following loop
     * would need to consider passing another check for
     * a separated big.matrix.  The following is required
     * for a separated big.matrix so it works in both scenarios
     */
     
    // get BigMatrix external pointer
    Rcpp::XPtr<BigMatrix> xpY = BigMatrixXPtr(Y);
    // create MatrixAccessor object to imput values
    MatrixAccessor<double> Y_BM(*xpY);
    
    // memcpy elements to the big.matrix object
    for (arma::uword i = 0; i < MM; i++){
      std::memcpy(Y_BM[i], Ym.colptr(i), NN*sizeof(double));
    }
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
//// [[Rcpp::export]]
//SEXP
//dgeqrf_wrapper (SEXP M, SEXP N, SEXP Y, SEXP LDA, SEXP TAU, SEXP WORK,
//                SEXP LWORK, SEXP INFO, SEXP A_isBM, SEXP TAU_isBM, 
//                SEXP WORK_isBM)

// [[Rcpp::export]]
SEXP
dgeqrf_wrapper (SEXP Y, SEXP Q, SEXP R)
{
//  #ifdef NON_R_BLAS
//    SEXP ans, Tr;
//    double *pY;
//    double *pTAU = make_double_ptr (TAU, TAU_isBM);
//    double *pWORK = make_double_ptr (WORK, WORK_isBM);
//    INT NN = (INT) * (DOUBLE_DATA (N));
//    INT MM = (INT) * (DOUBLE_DATA (M));
//    INT LDAi = (INT) * (DOUBLE_DATA (LDA));
//    INT LWORKi = (INT) * (DOUBLE_DATA (LWORK));
//    INT INFOi = (INT) * (DOUBLE_DATA (INFO));
//    
//    PROTECT(ans = Y);
//    PROTECT(Tr = Rf_allocVector(LGLSXP, 1));
//    LOGICAL(Tr)[0] = 1;
//    pY = make_double_ptr (Y, Tr);
//    
//    /* C-blas and REFBLAS untested */
//    /* An example of an alternate C-blas interface (e.g., ACML) */
//    #ifdef CBLAS
//      dgeqrf (MM, NN, pY, LDAi, pTAU, pWORK, LWORKi, INFOi);
//    #elif REFBLAS
//    /* Standard Fortran interface without underscoring */
//      int8_dgeqrf (&MM, &NN, pY, &LDAi, pTAU, pWORK, &LWORKi, &INFOi);
//    #else
//      #error "BLAS format not supported"
//    #endif
//    
//    Rf_unprotect(2);
//    return ans;
//  #else
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
    
    arma::uword QMM = Qm.n_cols;
    arma::uword QNN = Qm.n_rows;
    arma::uword RMM = Rm.n_cols;
    arma::uword RNN = Rm.n_rows;
    
     
    // get BigMatrix external pointer
    Rcpp::XPtr<BigMatrix> xpQ = BigMatrixXPtr(Q);
    Rcpp::XPtr<BigMatrix> xpR = BigMatrixXPtr(R);
    // create MatrixAccessor object to imput values
    MatrixAccessor<double> Q_BM(*xpQ);
    MatrixAccessor<double> R_BM(*xpR);
    
    // memcpy elements to the big.matrix object
    for (arma::uword i = 0; i < QMM; i++){
      std::memcpy(Q_BM[i], Qm.colptr(i), QNN*sizeof(double));
    }
    
    // memcpy elements to the big.matrix object
    for (arma::uword i = 0; i < RMM; i++){
      std::memcpy(R_BM[i], Rm.colptr(i), RNN*sizeof(double));
    }
    
    // return list of Q and R
    return Rcpp::List::create(Named("Q") = Q, Named("R") = R);
//  #endif
}


/* No specific BLAS functions for the following so all
 * are strictly armadillo.
 */

// element-wise matrix multiplcation
// [[Rcpp::export]]
SEXP
dgeemm_wrapper (SEXP X, SEXP Y, SEXP X_isBM){
  
    // convert to arma matrices
    const arma::mat Xm = ConvertToArma(X, X_isBM);
    // Y will always be a big.matrix
    arma::mat Ym = ConvertBMtoArma(Y);
      
    arma::uword MM = Xm.n_cols;
    arma::uword NN = Xm.n_rows;
    
    // RcppArmadillo Element-wise Matrix Multiplication
    Ym = Xm % Ym;
    
    // get BigMatrix external pointer
    Rcpp::XPtr<BigMatrix> xpY = BigMatrixXPtr(Y);
    // create MatrixAccessor object to imput values
    MatrixAccessor<double> Y_BM(*xpY);
    
    // memcpy elements to the big.matrix object
    for (arma::uword i = 0; i < MM; i++){
      std::memcpy(Y_BM[i], Ym.colptr(i), NN*sizeof(double));
    }
    return Y;
}


// element-wise matrix division
// [[Rcpp::export]]
SEXP
dgeemd_wrapper (SEXP X, SEXP Y, SEXP X_isBM){
    // convert to arma matrices
    const arma::mat Xm = ConvertToArma(X, X_isBM);
    // Y will always be a big.matrix
    arma::mat Ym = ConvertBMtoArma(Y);
      
    arma::uword MM = Xm.n_cols;
    arma::uword NN = Xm.n_rows;
    
    // RcppArmadillo Element-wise Matrix Division
    Ym = Xm / Ym;
    
    // get BigMatrix external pointer
    Rcpp::XPtr<BigMatrix> xpY = BigMatrixXPtr(Y);
    // create MatrixAccessor object to imput values
    MatrixAccessor<double> Y_BM(*xpY);
    
    // memcpy elements to the big.matrix object
    for (arma::uword i = 0; i < MM; i++){
      std::memcpy(Y_BM[i], Ym.colptr(i), NN*sizeof(double));
    }
    return Y;
}


// Scalar-matrix division
// [[Rcpp::export]]
SEXP
dgesmd_wrapper (SEXP A, SEXP Y, SEXP Y_isBM, int ALPHA_LHS) {

    INT ALPHA = (INT) * (DOUBLE_DATA(A));
  
    // Y will always be a big.matrix
    arma::mat Ym = ConvertBMtoArma(Y);
      
    arma::uword MM = Ym.n_cols;
    arma::uword NN = Ym.n_rows;
    
    // RcppArmadillo Element-wise Scalar Division
    if(ALPHA_LHS){
        Ym = ALPHA / Ym;
    }else{      
        Ym = Ym / ALPHA;
    }
    
    // get BigMatrix external pointer
    Rcpp::XPtr<BigMatrix> xpY = BigMatrixXPtr(Y);
    // create MatrixAccessor object to imput values
    MatrixAccessor<double> Y_BM(*xpY);
    
    // memcpy elements to the big.matrix object
    for (arma::uword i = 0; i < MM; i++){
      std::memcpy(Y_BM[i], Ym.colptr(i), NN*sizeof(double));
    }
    return Y;
}


// Power
// [[Rcpp::export]]
SEXP
dgepow_wrapper(SEXP EXP, SEXP Y) {
    // Y will always be a big.matrix
    arma::mat Ym = ConvertBMtoArma(Y);
      
    arma::uword MM = Ym.n_cols;
    arma::uword NN = Ym.n_rows;
    
    // RcppArmadillo Element-wise Scalar Division
    Ym = pow(Ym, Rcpp::as<double>(EXP));
    
    // get BigMatrix external pointer
    Rcpp::XPtr<BigMatrix> xpY = BigMatrixXPtr(Y);
    // create MatrixAccessor object to imput values
    MatrixAccessor<double> Y_BM(*xpY);
    
    // memcpy elements to the big.matrix object
    for (arma::uword i = 0; i < MM; i++){
      std::memcpy(Y_BM[i], Ym.colptr(i), NN*sizeof(double));
    }
    return Y;
}


// common logarithm
// [[Rcpp::export]]
SEXP
dgeclog_wrapper(SEXP Y) {
    // Y will always be a big.matrix
    arma::mat Ym = ConvertBMtoArma(Y);
      
    arma::uword MM = Ym.n_cols;
    arma::uword NN = Ym.n_rows;
    
    // RcppArmadillo Element-wise Scalar Division
    Ym = log10(Ym);
    
    // get BigMatrix external pointer
    Rcpp::XPtr<BigMatrix> xpY = BigMatrixXPtr(Y);
    // create MatrixAccessor object to imput values
    MatrixAccessor<double> Y_BM(*xpY);
    
    // memcpy elements to the big.matrix object
    for (arma::uword i = 0; i < MM; i++){
      std::memcpy(Y_BM[i], Ym.colptr(i), NN*sizeof(double));
    }
    return Y;
}

// base logarithm
// [[Rcpp::export]]
SEXP
dgelog_wrapper(SEXP BASE, SEXP Y) {
    // Y will always be a big.matrix
    arma::mat Ym = ConvertBMtoArma(Y);
      
    arma::uword MM = Ym.n_cols;
    arma::uword NN = Ym.n_rows;
    
    // RcppArmadillo Element-wise Scalar Division
    Ym = log10(Ym)/log10(Rcpp::as<double>(BASE));
    
    // get BigMatrix external pointer
    Rcpp::XPtr<BigMatrix> xpY = BigMatrixXPtr(Y);
    // create MatrixAccessor object to imput values
    MatrixAccessor<double> Y_BM(*xpY);
    
    // memcpy elements to the big.matrix object
    for (arma::uword i = 0; i < MM; i++){
      std::memcpy(Y_BM[i], Ym.colptr(i), NN*sizeof(double));
    }
    return Y;
}

// Exponential function
// [[Rcpp::export]]
SEXP
dgeexp_wrapper(SEXP Y) {
    // Y will always be a big.matrix
    arma::mat Ym = ConvertBMtoArma(Y);
      
    arma::uword MM = Ym.n_cols;
    arma::uword NN = Ym.n_rows;
    
    // RcppArmadillo Element-wise Scalar Division
    Ym = exp(Ym);
    
    // get BigMatrix external pointer
    Rcpp::XPtr<BigMatrix> xpY = BigMatrixXPtr(Y);
    // create MatrixAccessor object to imput values
    MatrixAccessor<double> Y_BM(*xpY);
    
    // memcpy elements to the big.matrix object
    for (arma::uword i = 0; i < MM; i++){
      std::memcpy(Y_BM[i], Ym.colptr(i), NN*sizeof(double));
    }
    return Y;
}

// hyperbolic tangent
// [[Rcpp::export]]
SEXP
dgetanh_wrapper(SEXP Y) {
    // Y will always be a big.matrix
    arma::mat Ym = ConvertBMtoArma(Y);
      
    arma::uword MM = Ym.n_cols;
    arma::uword NN = Ym.n_rows;
    
    // RcppArmadillo Element-wise Scalar Division
    Ym = tanh(Ym);
    
    // get BigMatrix external pointer
    Rcpp::XPtr<BigMatrix> xpY = BigMatrixXPtr(Y);
    // create MatrixAccessor object to imput values
    MatrixAccessor<double> Y_BM(*xpY);
    
    // memcpy elements to the big.matrix object
    for (arma::uword i = 0; i < MM; i++){
      std::memcpy(Y_BM[i], Ym.colptr(i), NN*sizeof(double));
    }
    return Y;
}

// hyperbolic cosine
// [[Rcpp::export]]
SEXP
dgecosh_wrapper(SEXP Y) {
    // Y will always be a big.matrix
    arma::mat Ym = ConvertBMtoArma(Y);
      
    arma::uword MM = Ym.n_cols;
    arma::uword NN = Ym.n_rows;
    
    // RcppArmadillo Element-wise Scalar Division
    Ym = cosh(Ym);
    
    // get BigMatrix external pointer
    Rcpp::XPtr<BigMatrix> xpY = BigMatrixXPtr(Y);
    // create MatrixAccessor object to imput values
    MatrixAccessor<double> Y_BM(*xpY);
    
    // memcpy elements to the big.matrix object
    for (arma::uword i = 0; i < MM; i++){
      std::memcpy(Y_BM[i], Ym.colptr(i), NN*sizeof(double));
    }
    return Y;
}

// hyperbolic sine
// [[Rcpp::export]]
SEXP
dgesinh_wrapper(SEXP Y) {
    // Y will always be a big.matrix
    arma::mat Ym = ConvertBMtoArma(Y);
      
    arma::uword MM = Ym.n_cols;
    arma::uword NN = Ym.n_rows;
    
    // RcppArmadillo Element-wise Scalar Division
    Ym = sinh(Ym);
    
    // get BigMatrix external pointer
    Rcpp::XPtr<BigMatrix> xpY = BigMatrixXPtr(Y);
    // create MatrixAccessor object to imput values
    MatrixAccessor<double> Y_BM(*xpY);
    
    // memcpy elements to the big.matrix object
    for (arma::uword i = 0; i < MM; i++){
      std::memcpy(Y_BM[i], Ym.colptr(i), NN*sizeof(double));
    }
    return Y;
}

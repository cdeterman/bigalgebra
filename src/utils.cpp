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

#include "bigmatrix_arma.h"

using namespace Rcpp;

// [[Rcpp::export]]
SEXP
all_equal_cpp (SEXP X_, SEXP Y_, SEXP X_isBM, SEXP Y_isBM, SEXP tol_)
{
  arma::mat X = ConvertToArma(X_, X_isBM);
  arma::mat Y = ConvertToArma(Y_, Y_isBM);
  double tol = as<double>(tol_);
  double check = arma::sum(arma::sum(X-Y));
  
  if(check > tol){
      return wrap(false);
  }else{
      return wrap(true);
  }
}

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

#include "bigmatrix_arma.h"

// [[Rcpp::export]]
SEXP
all_equal_cpp (SEXP X_, SEXP Y_, SEXP tol_)
{
  arma::mat X = ConvertBMtoArma(X_);
  arma::mat Y = ConvertBMtoArma(Y_);
  double tol = Rcpp::as<double>(tol_);
  double check = arma::sum(arma::sum(X-Y));
  
  if(check > tol){
      return Rcpp::wrap(false);
  }else{
      return Rcpp::wrap(true);
  }
}
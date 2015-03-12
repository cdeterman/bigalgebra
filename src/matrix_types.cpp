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
cpp_isDiagonal(SEXP Y, SEXP Y_isBM){
  
  const arma::mat Ym = ConvertToArma(Y, Y_isBM);
  
  // lower matrix
  const arma::vec Yl = vectorise(arma::trimatl(Ym));
  // upper matrix
  const arma::vec Yu = vectorise(arma::trimatu(Ym));

  bool out;
  out = all(Yl == Yu);
  return wrap(out);
}


// [[Rcpp::export]]
LogicalVector
cpp_isTriangular(SEXP Y, SEXP Y_isBM, LogicalVector upper_){
  
  //bool upper = as<bool>(upper_);
  LogicalVector out = LogicalVector::create(true);
  const arma::mat Ym = ConvertToArma(Y, Y_isBM);
  const arma::uword nc = Ym.n_cols;
  
  arma::mat Z(nc,nc,arma::fill::zeros);
  arma::mat X(nc,nc,arma::fill::zeros);
  
  arma::vec idx = arma::linspace<arma::mat>(1,nc,nc);
  
  X.each_col() += idx;
  Z.each_row() += trans(idx);
  
  if(any(is_na(upper_))){
      const arma::vec u = Ym.elem(arma::find(Z<X));
      const arma::vec l = Ym.elem(arma::find(Z>X));
      bool upper = all(u==0);
      bool lower = all(l==0);
      
      if(upper && lower){
        return out;
      }else if(upper){
        out.attr("kind") = "U";
        return out;
      }else if(lower){
        out.attr("kind") = "L";
        return out;
      }else{
        out[0] = false;
        return out;
      }
  }else{
    bool upper = as<bool>(upper_);
    if(!upper){
      const arma::vec l = Ym.elem(arma::find(Z>X));
      if(all(l==0)){
        out.attr("kind") = "L";
        return out;
      }else{
        out[0] = false;
        return out;
      }
    }else if(upper){
        const arma::vec u = Ym.elem(arma::find(Z<X));
        if(all(u==0)){
          out.attr("kind") = "U";
          return out;
        }else{
          out[0] = false;
          return out;
        }
    }else{
      stop("Unrecognized upper selection");
    }
  }
  return LogicalVector::create(false);
}

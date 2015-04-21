#ifndef BIGMATRIX_ARMA
#define BIGMATRIX_ARMA

// simple function to get access to big.matrix objects from R
inline 
Rcpp::XPtr<BigMatrix> BigMatrixXPtr(SEXP A){
    // declare as S4 object
    Rcpp::S4 As4(A);
    // pull address slot
    SEXP A_address = As4.slot("address");
    // declare as external pointer
    Rcpp::XPtr<BigMatrix> xpA(A_address);
    return xpA;
 }
 
 /* Directly convert big.matrix to arma object
 * NOTE - this shares the memory so anything
 * done to the arma object is automatically changed
 * in the big.matrix object.
 */
inline
arma::mat ConvertBMtoArma(SEXP A)
{
  Rcpp::XPtr<BigMatrix> xpA = BigMatrixXPtr(A);
        
  arma::mat Am = arma::mat( (double*) xpA->matrix(),
              xpA->nrow(),
              xpA->ncol(),
              false);
  return Am;
}


#endif
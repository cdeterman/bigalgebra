
#' @import assertive
assert_is_bigmatrix <- function (x) 
{
  assertive:::assert_engine(x, is.big.matrix, "x is not a big.matrix")
}

#' @export
assert_is_positive_definite <- function(x)
{
  assertive:::assert_engine(x, isPositiveDefinite, 
                            "matrix is not positive-definite")
}
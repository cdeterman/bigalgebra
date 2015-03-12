
#' @import assertive
assert_is_bigmatrix <- function (x) 
{
  assertive:::assert_engine(x, is.big.matrix, "x is not a big.matrix")
}

assert_is_square <- function(x)
{
  assertive:::assert_engine(x, isSquare, 
                            "matrix is not square")
}

#' @export
assert_is_positive_definite <- function(x)
{
  assertive:::assert_engine(x, isPositiveDefinite, 
                            "matrix is not positive-definite")
}

#' @export
assert_is_positive_semi_definite <- function(x)
{
  assertive:::assert_engine(x, isPositiveSemiDefinite, 
                            "matrix is not positive semi-definite")
}

#' @export
assert_is_negative_definite <- function(x)
{
  assertive:::assert_engine(x, isNegativeDefinite, 
                            "matrix is not negative-definite")
}

#' @export
assert_is_negative_semi_definite <- function(x)
{
  assertive:::assert_engine(x, isNegativeSemiDefinite, 
                            "matrix is not negative semi-definite")
}


#' @import assertive
assert_is_bigmatrix <- function (x) 
{
  if(!is.big.matrix(x)){
    stop("x is not a big.matrix")
  }
}

assert_is_square <- function(x)
{
  if(!isSquare(x)){
    stop("matrix is not square")
  }
}

#' @export
assert_is_positive_definite <- function(x)
{
  if(!isPositiveDefinite(x)){
    stop("matrix is not positive-definite")
  }
}

#' @export
assert_is_positive_semi_definite <- function(x)
{
  if(!isPositiveSemiDefinite(x)){
    stop("matrix is not positive semi-definite")
  }
}

#' @export
assert_is_negative_definite <- function(x)
{
  if(!isNegativeDefinite(x)){
    stop("matrix is not negative-definite")
  }
}

#' @export
assert_is_negative_semi_definite <- function(x)
{
  if(!isNegativeSemiDefinite(x)){
    stop("matrix is not negative semi-definite")
  }
}

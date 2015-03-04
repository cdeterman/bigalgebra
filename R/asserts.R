
#' @import assertive
assert_is_bigmatrix <- function (x) 
{
  assertive:::assert_engine(x, is.big.matrix)
}

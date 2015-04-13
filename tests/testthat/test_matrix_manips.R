library("bigalgebra")
context("Matrix Manipulation Tests")

mat <- matrix(1:12, ncol = 4)
mat2 <- matrix(rnorm(12), ncol = 4)
mat3 <- matrix(rnorm(15), ncol = 5)

storage.mode(mat) <-"double"
storage.mode(mat2) <- "double"
storage.mode(mat3) <- "double"

bm <- as.big.matrix(mat, type="double")
bm2 <- as.big.matrix(mat2, type="double")
bm3 <- as.big.matrix(mat3, type="double")

# Diagonal Matrix
dm <- diag(rep(1,3))
dBM <- as.big.matrix(dm, type="double")

# triangular matrices
utm <- matrix(seq(25), 5); utm[!upper.tri(utm)] <- 0
diag(utm) <- 1
utmBM <- as.big.matrix(utm, type="double")
ltm <- matrix(seq(25), 5); ltm[!lower.tri(ltm)] <- 0
ltmBM <- as.big.matrix(ltm, type="double")

test_that("Transposition successful",{
  expect_true(class(t(bm)) == "big.matrix")
  expect_equivalent(t(mat), t(bm)[])
})

test_that("all.equal method successful",{
  expect_true(all.equal(bm, bm))
  expect_false(all.equal(bm, bm2))
  expect_false(all.equal(bm, bm3))
})

test_that("isDiagonal works", {  
  expect_true(isDiagonal(dm))
  expect_true(isDiagonal(dBM))
  expect_false(isDiagonal(mat))
  expect_false(isDiagonal(bm))
  expect_false(isDiagonal(utm))
  expect_false(isDiagonal(utmBM))
  expect_false(isDiagonal(ltm))
  expect_false(isDiagonal(ltmBM))
})

test_that("isTriangular works", {
#   print(dm)
  expect_true(isTriangular(dm))
  expect_true(isTriangular(dBM))
  expect_true(isTriangular(utm, upper=TRUE), 
              "upper triangular matrix should be true")
  expect_true(isTriangular(utmBM, upper=TRUE), 
              "upper triangular big.matrix should be true")
  expect_true(isTriangular(ltm, upper=FALSE), 
              "lower triangular matrix should be true")
  expect_true(isTriangular(ltmBM, upper=FALSE), 
              "lower triangular big.matrix should be true")
  expect_true(isTriangular(utm, upper=NA), 
              "upper triangular matrix should be true")
  expect_true(isTriangular(utmBM, upper=NA), 
              "upper triangular big.matrix should be true")
  expect_true(isTriangular(ltm, upper=NA), 
              "lower triangular matrix should be true")
  expect_true(isTriangular(ltmBM, upper=NA), 
              "lower triangular big.matrix should be true")
  expect_false(isTriangular(utm, upper=FALSE), 
               "upper triangular matrix should be false")
  expect_false(isTriangular(utmBM, upper=FALSE), 
               "upper triangular big.matrix should be false")
  expect_false(isTriangular(ltm, upper=TRUE), 
               "lower triangular matrix should be false")
  expect_false(isTriangular(ltmBM, upper=TRUE), 
               "lower triangular big.matrix should be false")
  expect_true(attr(isTriangular(utm, TRUE), "kind") == "U", 
              "incorrect attribute assignment")
  expect_true(attr(isTriangular(utmBM, TRUE), "kind") == "U", 
              "incorrect big.matrix attribute assignment")
  expect_true(attr(isTriangular(ltm, FALSE), "kind") == "L", 
              "incorrect attribute assignment")
  expect_true(attr(isTriangular(ltmBM, FALSE), "kind") == "L", 
              "incorrect big.matrix attribute assignment")
  expect_true(attr(isTriangular(utm, NA), "kind") == "U", 
              "incorrect attribute assignment")
  expect_true(attr(isTriangular(utmBM, NA), "kind") == "U", 
              "incorrect big.matrix attribute assignment")
  expect_true(attr(isTriangular(ltm, NA), "kind") == "L", 
              "incorrect attribute assignment") 
  expect_true(attr(isTriangular(ltmBM, NA), "kind") == "L", 
              "incorrect big.matrix attribute assignment")  
})

test_that("isSymmetric method works",{
  # create symmetric matrix
  s <- matrix(seq(25), 5)
  s[lower.tri(s)] = t(s)[lower.tri(s)]
  symBM <- as.big.matrix(s, type="double")
  
  expect_true(isSymmetric(s))
  expect_true(isSymmetric(symBM))
})

test_that("isPositiveDefinite & isPositiveSemiDefinite methods work",{
  # create positive definite matrix
  pd <- matrix( c( 2, -1, 0, -1 , 2, -1, 0, -1, 2 ), nrow=3, byrow=TRUE)
  # positive semi-definite matrix
  psd <- matrix( c( 2, -1, 2, -1, 2, -1, 2, -1, 2 ), nrow=3, byrow=TRUE )

  # matrixcalc::is.positive.definite(pd)
  pdBM <- as.big.matrix(pd,type="double")
  psdBM <- as.big.matrix(psd,type="double")
  
  expect_true(isPositiveDefinite(pd))
  expect_true(isPositiveDefinite(pdBM))
  expect_error(isPositiveDefinite(mat))
  expect_error(isPositiveDefinite(bm))
  expect_false(isPositiveDefinite(psd))
  expect_false(isPositiveDefinite(psdBM))
  expect_true(isPositiveSemiDefinite(psd))
  expect_true(isPositiveSemiDefinite(psdBM))
})

test_that("isNegativeDefinite & isNegativeSemiDefinite methods works",{
  # create negative definite matrix
  nd <- matrix( c( -2, 1, 0, 1, -2, 1, 0, 1, -2 ), nrow=3, byrow=TRUE )
  # negative semi-definite matrix
  nsd <- diag(c(-5,0,-1))

  # matrixcalc::is.positive.definite(nd)
  ndBM <- as.big.matrix(nd,type="double")
  nsdBM <- as.big.matrix(nsd,type="double")
  
  expect_true(isNegativeDefinite(nd))
  expect_true(isNegativeDefinite(ndBM))
  expect_error(isNegativeDefinite(mat))
  expect_error(isNegativeDefinite(bm))
  expect_false(isNegativeDefinite(nsd))
  expect_false(isNegativeDefinite(nsdBM))
  expect_true(isNegativeSemiDefinite(nsd))
  expect_true(isNegativeSemiDefinite(nsdBM))
})


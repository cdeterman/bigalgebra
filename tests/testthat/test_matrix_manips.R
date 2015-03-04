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

test_that("Transposition successful",{
  expect_true(class(t(bm)) == "big.matrix")
  expect_equivalent(t(mat), t(bm)[])
})

test_that("all.equal method successful",{
  expect_true(all.equal(bm, bm))
  expect_false(all.equal(bm, bm2))
  expect_false(all.equal(bm, bm3))
})

test_that("isSymmetric method works",{
  # create symmetric matrix
  s <- matrix(seq(25), 5)
  s[lower.tri(s)] = t(s)[lower.tri(s)]
  symBM <- as.big.matrix(s, type="double")
  
  expect_true(isSymmetric(s))
  expect_true(isSymmetric(symBM))
})

test_that("isPositiveDefinite method works",{
  # create positive definite matrix
  pd <- matrix(c(2,-1,0,-1,2,-1,0,-1,2), nrow=3, byrow=TRUE)
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
})


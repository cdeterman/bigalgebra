library("bigalgebra")
context("Linear Algebra Tests")

mat <- matrix(1:9, 
              ncol = 3)
mat2 <- matrix(1:9, 
              ncol = 3)
mat3 <- matrix(1:12, 
               ncol = 4)


# set as type double
storage.mode(mat) <- "double"
storage.mode(mat2) <- "double"
storage.mode(mat3) <- "double"

bm <- as.big.matrix(mat, type="double")
bm2 <- as.big.matrix(mat2, type="double")
bm3 <- as.big.matrix(mat3, type="double")

test_that("matrix multiplication successful",{
  R_mat <- mat %*% mat2
  BM_mat <- bm %*% bm2
  R_BM_mat <- bm %*% mat2
  BM_R_mat <- mat %*% bm2
  expect_equivalent(R_mat, BM_mat[,])
  expect_equivalent(R_mat, R_BM_mat[,])
  expect_equivalent(R_mat, BM_R_mat[,])
  expect_error(bm3 %*% bm2)
})

test_that("matrix addition successful", {
  R_mat_s <- mat + 1
  R_mat_m <- mat + mat
  BM_mat_s <- bm + 1
  BM_mat_m <- bm + bm
  expect_equivalent(R_mat_s, BM_mat_s[,])
  expect_equivalent(R_mat_m, BM_mat_m[,])
})

test_that("matrix subtraction successful", {
  
  R_mat_m <- mat - mat
  BM_mat_m <- bm - bm
  expect_equivalent(R_mat_m, BM_mat_m[,])
})

test_that("matrix scalar subtraction successful", {
  R_mat_s <- mat - 1
  R_mat_s2 <- 1 - mat
  BM_mat_s <- bm - 1
  BM_mat_s2 <- 1 - bm
  
  expect_equivalent(R_mat_s, BM_mat_s[,])
  expect_equivalent(R_mat_s2, BM_mat_s2[,])
})

test_that("matrix unary subtraction successful", {
  R_mat_s <- -mat
  BM_mat_s <- -bm
  
  expect_equivalent(R_mat_s, BM_mat_s[,])
})

test_that("scalar matrix mutliplication successful", {
  R_mat <- 3 * mat
  BM_mat_sm <- 3 * bm
  BM_mat_ms <- bm * 3
  expect_equivalent(R_mat, BM_mat_sm[,])
  expect_equivalent(R_mat, BM_mat_ms[,])
})

test_that("matrix element-wise multiplication successful", {
  R_mat <- mat * mat2
  BM_mat <- bm * bm2
  R_BM_mat <- bm * mat2
  BM_R_mat <- mat * bm2
  expect_equivalent(R_mat, BM_mat[,])
  expect_equivalent(R_mat, R_BM_mat[,])
  expect_equivalent(R_mat, BM_R_mat[,])
  expect_error(bm3 * bm2)
})

test_that("matrix element-wise division successful", {
  R_mat <- mat / mat2
  R_mat_s_num <- 1 / mat
  R_mat_s_denom <- mat / 2
  BM_mat <- bm / bm2
  BM_mat_s_num <- 1 / bm
  BM_mat_s_denom <- bm / 2
  R_BM_mat <- bm / mat2
  BM_R_mat <- mat / bm2
  expect_equivalent(R_mat, BM_mat[,])
  expect_equivalent(R_mat, R_BM_mat[,])
  expect_equivalent(R_mat, BM_R_mat[,])
  expect_equivalent(R_mat_s_num, BM_mat_s_num[,])
  expect_equivalent(R_mat_s_denom, BM_mat_s_denom[,])
  expect_error(bm3 / bm2)
})

test_that("matrix power is successful", {
  R_mat <- mat ^ 2
  BM_mat <- bm ^ 2
  expect_equivalent(R_mat, BM_mat[,])
})

test_that("math functions successful", {
  R_mat_sqrt <- sqrt(mat)
  R_mat_exp <- exp(mat)
  R_mat_log10 <- log10(mat)
  R_mat_sinh <- sinh(mat)
  R_mat_cosh <- cosh(mat)
  R_mat_tanh <- tanh(mat)
  BM_mat_sqrt <- sqrt(bm)
  BM_mat_exp <- exp(bm)
  BM_mat_log10 <- log10(bm)
  BM_mat_sinh <- sinh(bm)
  BM_mat_cosh <- cosh(bm)
  BM_mat_tanh <- tanh(bm)
  
  expect_equivalent(R_mat_sqrt, BM_mat_sqrt[,])
  expect_equivalent(R_mat_exp, BM_mat_exp[,])
  expect_equivalent(R_mat_log10, BM_mat_log10[,])
  expect_equivalent(R_mat_sinh, BM_mat_sinh[,])
  expect_equivalent(R_mat_cosh, BM_mat_cosh[,])
  expect_equivalent(R_mat_tanh, BM_mat_tanh[,])
})

test_that("QR matrix decomposition successful", {
  BM_QR <- qr(bm)
  expect_true(all(sapply(BM_QR, function(x) class(x) == "big.matrix")))
  expect_equivalent(BM_QR$Q[] %*% BM_QR$R[], bm[])
})

test_that("Choleski matrix decomposition successful", {
  m_chol <- matrix(c(2,-1,0,-1,2,-1,0,-1,2), nrow=3, byrow=TRUE)
  bm_chol <- as.big.matrix(m_chol)
  M_CHOL <- chol(m_chol)
  BM_CHOL <- chol(bm_chol)
  expect_true(class(BM_CHOL) == "big.matrix")
  expect_equivalent(M_CHOL, BM_CHOL[])
})

test_that("eigen method works", {
  # create symmetric matrix
  s <- matrix(seq(25), 5)
  s[lower.tri(s)] = t(s)[lower.tri(s)]
  symBM <- as.big.matrix(s, type="double")
  
  mat_eig <- eigen(s, only.values=TRUE)$values
  BM_eig <- eigen(symBM, only.values=TRUE)$values
  mat_eig_full <- eigen(s, only.values=FALSE)
  BM_eig_full <- eigen(symBM, only.values=FALSE)
  
  expect_equivalent(mat_eig, BM_eig)
  expect_equivalent(mat_eig_full$values, BM_eig_full$values)
  expect_true(class(BM_eig_full$vectors) == "big.matrix")
  # eigen vectors signs are irrelevant so ignore with abs
  expect_equivalent(abs(mat_eig_full$vectors), abs(BM_eig_full$vectors[]))
})

test_that("crossprod method", {
  
  X <- matrix(rnorm(10), nrow=2)
  Y <- matrix(rnorm(10), nrow=2)
  Z <- matrix(rnorm(10), nrow=5)
  
  C <- crossprod(X,Y)
  Cs <- crossprod(X)
  
  Xbm <- as.big.matrix(X)
  Ybm <- as.big.matrix(Y)
  Zbm <- as.big.matrix(Z)
  
  Cbm <- crossprod(Xbm, Ybm)
  Csbm <- crossprod(Xbm)
  
  expect_is(Cbm, "big.matrix")
  expect_equal(Cbm[,], C, tolerance=.Machine$double.eps ^ 0.5, 
               info="big.matrix elements not equivalent")  
  expect_equal(Csbm[,], Cs, tolerance=.Machine$double.eps ^ 0.5, 
               info="big.matrix elements not equivalent") 
  expect_error(crossprod(Xbm, Zbm))
})

test_that("tcrossprod method", {
  
  X <- matrix(rnorm(10), nrow=5)
  Y <- matrix(rnorm(10), nrow=5)
  Z <- matrix(rnorm(15), nrow=5)
  
  C <- tcrossprod(X,Y)
  Cs <- tcrossprod(X)
  
  Xbm <- as.big.matrix(X)
  Ybm <- as.big.matrix(Y)
  Zbm <- as.big.matrix(Z)
  
  Cbm <- tcrossprod(Xbm, Ybm)
  Csbm <- tcrossprod(Xbm)
  
  expect_is(Cbm, "big.matrix")
  expect_equal(Cbm[,], C, tolerance=.Machine$double.eps ^ 0.5, 
               info="big.matrix elements not equivalent")  
  expect_equal(Csbm[,], Cs, tolerance=.Machine$double.eps ^ 0.5, 
               info="big.matrix elements not equivalent") 
  expect_error(tcrossprod(Xbm, Zbm))
})



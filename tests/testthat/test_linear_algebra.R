library("bigalgebra")
context("linear algebra tests")

mat <- matrix(1:9, 
              ncol = 3, 
              nrow = 3, 
              dimnames = list(letters[1:3], LETTERS[1:3]))
mat2 <- matrix(1:9, 
              ncol = 3, 
              nrow = 3, 
              dimnames = list(letters[1:3], LETTERS[1:3]))

# set as type double
storage.mode(mat) <- "double"
storage.mode(mat2) <- "double"

bm <- as.big.matrix(mat, type="double")
bm2 <- as.big.matrix(mat2, type="double")

C <- bigalgebra:::anon_matrix(3, 3, "double")

bigalgebra:::dgemm_wrapper2(A=bm@address, 
                            B=bm2@address, 
                            C=C@address,
                            TRUE, TRUE, TRUE)[,]
# bm[,]%*%bm2[,]

test_that("matrix multiplication successful",{
  R_mat <- mat %*% mat2
  BM_mat <- bm %*% bm2
  R_BM_mat <- bm %*% mat2
  BM_R_mat <- mat %*% bm2
  expect_equivalent(R_mat, BM_mat[,])
  expect_equivalent(R_mat, R_BM_mat[,])
  expect_equivalent(R_mat, BM_R_mat[,])
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
  R_mat_s <- mat - 1
  R_mat_m <- mat - mat
  BM_mat_s <- bm - 1
  BM_mat_m <- bm - bm
  expect_equivalent(R_mat_s, BM_mat_s[,])
  expect_equivalent(R_mat_m, BM_mat_m[,])
})

test_that("matrix element-wise multiplication successful", {
  R_mat <- mat * mat2
  BM_mat <- bm * bm2
  R_BM_mat <- bm * mat2
  BM_R_mat <- mat * bm2
  expect_equivalent(R_mat, BM_mat[,])
  expect_equivalent(R_mat, R_BM_mat[,])
  expect_equivalent(R_mat, BM_R_mat[,])
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
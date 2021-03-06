# This file was generated by Rcpp::compileAttributes
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @useDynLib bigalgebra
#' @importFrom Rcpp evalCpp
NULL

dgemm_wrapper <- function(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC, A_isBM, B_isBM, C_isBM, C_offset) {
    .Call('bigalgebra_dgemm_wrapper', PACKAGE = 'bigalgebra', TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC, A_isBM, B_isBM, C_isBM, C_offset)
}

daxpy_wrapper <- function(N, A, X, Y, X_isBM) {
    .Call('bigalgebra_daxpy_wrapper', PACKAGE = 'bigalgebra', N, A, X, Y, X_isBM)
}

daxpy_unary_wrapper <- function(X, X_isBM) {
    invisible(.Call('bigalgebra_daxpy_unary_wrapper', PACKAGE = 'bigalgebra', X, X_isBM))
}

dpotrf_wrapper <- function(UPLO, N, A, LDA, INFO, A_isBM) {
    .Call('bigalgebra_dpotrf_wrapper', PACKAGE = 'bigalgebra', UPLO, N, A, LDA, INFO, A_isBM)
}

dadd_wrapper <- function(ALPHA, Y, S, Y_isBM) {
    invisible(.Call('bigalgebra_dadd_wrapper', PACKAGE = 'bigalgebra', ALPHA, Y, S, Y_isBM))
}

dgeqrf_wrapper <- function(Y, Q, R) {
    .Call('bigalgebra_dgeqrf_wrapper', PACKAGE = 'bigalgebra', Y, Q, R)
}

dgeemm_wrapper <- function(X, Y, X_isBM) {
    .Call('bigalgebra_dgeemm_wrapper', PACKAGE = 'bigalgebra', X, Y, X_isBM)
}

dgeemd_wrapper <- function(X, Y, X_isBM) {
    .Call('bigalgebra_dgeemd_wrapper', PACKAGE = 'bigalgebra', X, Y, X_isBM)
}

dgesmd_wrapper <- function(A, Y, Y_isBM, ALPHA_LHS) {
    .Call('bigalgebra_dgesmd_wrapper', PACKAGE = 'bigalgebra', A, Y, Y_isBM, ALPHA_LHS)
}

#' @name t_wrapper
#' @title Big Matrix Transposition
#' @param X big.matrix object
#' @param Y big.matrix object (to be filled)
#' @param LOW_MEM boolean option to choose slower, low memory option
t_wrapper <- function(X, Y) {
    .Call('bigalgebra_t_wrapper', PACKAGE = 'bigalgebra', X, Y)
}

#' @name t_inplace_wrapper
#' @title Big Matrix In Place Transposition
#' @param X big.matrix object
#' @param LOW_MEM boolean option to choose slower, low memory option
t_inplace_wrapper <- function(X, LOW_MEM) {
    .Call('bigalgebra_t_inplace_wrapper', PACKAGE = 'bigalgebra', X, LOW_MEM)
}

cpp_bm_crossprod <- function(X, Y, Z, X_isBM, Y_isBM, Z_isBM) {
    invisible(.Call('bigalgebra_cpp_bm_crossprod', PACKAGE = 'bigalgebra', X, Y, Z, X_isBM, Y_isBM, Z_isBM))
}

cpp_bm_tcrossprod <- function(X, Y, Z, X_isBM, Y_isBM, Z_isBM) {
    invisible(.Call('bigalgebra_cpp_bm_tcrossprod', PACKAGE = 'bigalgebra', X, Y, Z, X_isBM, Y_isBM, Z_isBM))
}

dgepow_wrapper <- function(EXP, Y) {
    .Call('bigalgebra_dgepow_wrapper', PACKAGE = 'bigalgebra', EXP, Y)
}

dgeclog_wrapper <- function(Y) {
    .Call('bigalgebra_dgeclog_wrapper', PACKAGE = 'bigalgebra', Y)
}

dgelog_wrapper <- function(BASE, Y) {
    .Call('bigalgebra_dgelog_wrapper', PACKAGE = 'bigalgebra', BASE, Y)
}

dgeexp_wrapper <- function(Y) {
    .Call('bigalgebra_dgeexp_wrapper', PACKAGE = 'bigalgebra', Y)
}

dgetanh_wrapper <- function(Y) {
    .Call('bigalgebra_dgetanh_wrapper', PACKAGE = 'bigalgebra', Y)
}

dgecosh_wrapper <- function(Y) {
    .Call('bigalgebra_dgecosh_wrapper', PACKAGE = 'bigalgebra', Y)
}

dgesinh_wrapper <- function(Y) {
    .Call('bigalgebra_dgesinh_wrapper', PACKAGE = 'bigalgebra', Y)
}

eigen_wrapper <- function(X, EIG_VECS, only_values) {
    .Call('bigalgebra_eigen_wrapper', PACKAGE = 'bigalgebra', X, EIG_VECS, only_values)
}

cpp_isDiagonal <- function(Y, Y_isBM) {
    .Call('bigalgebra_cpp_isDiagonal', PACKAGE = 'bigalgebra', Y, Y_isBM)
}

cpp_isTriangular <- function(Y, Y_isBM, upper_) {
    .Call('bigalgebra_cpp_isTriangular', PACKAGE = 'bigalgebra', Y, Y_isBM, upper_)
}

all_equal_cpp <- function(X_, Y_, X_isBM, Y_isBM, tol_) {
    .Call('bigalgebra_all_equal_cpp', PACKAGE = 'bigalgebra', X_, Y_, X_isBM, Y_isBM, tol_)
}


# XXX TODO
#
# Each routine below may return either a big matrix or an R matrix
# depending on how they are called (or a set of such matrices for
# the decomposition methods). In each case the returned value is
# allocated by the routine.
#
# When big matrices are returned, we need to register a finalizer
# with the big matrix that removes its allocated files when the
# garbage collector is called on the big matrix, to clean up.
# Add this!
#


# Matrix Multiply
# C := ALPHA * op(A) * op(B) + BETA * C
# This is function provides dgemm functionality.
dgemm = function(TRANSA='N', TRANSB='N', M=NULL, N=NULL, K=NULL,
  ALPHA=1, A, LDA=NULL, B, LDB=NULL, BETA=0, C, LDC=NULL, COFF=0) 
{
  
  if(ncol(A) != nrow(B)){
    stop("matrices are not conformable") 
  }
  
  A.is.bm = check_matrix(A)
  B.is.bm = check_matrix(B)
  
# The matrices look OK.  Now, if they haven't been specified, let's
# specify some reasonable dimension information.
  if ( is.null(M) )
  {
    M = ifelse ( is_transposed(TRANSA), ncol(A), nrow(A) )
  }
  if ( is.null(N) ) 
  {
    N = ifelse ( is_transposed(TRANSB), nrow(B), ncol(B) )
  }
  if ( is.null(K) )
  {
    K = ifelse ( is_transposed(TRANSA), nrow(A), ncol(A) )
  }
  if ( is.null(LDA) ) LDA = ifelse (is_transposed(TRANSA), K, M)
  if ( is.null(LDB) ) LDB = ifelse (is_transposed(TRANSB), N, K) 
  if ( is.null(LDC) ) LDC = M

  # Default to big matrix output
  if(missing(C)) C = anon_matrix(M, N)
  C.is.bm = "big.matrix" %in% class(C)

  dgemm_wrapper(as.character(TRANSA), as.character(TRANSB),
    as.double(M), as.double(N), as.double(K), as.double(ALPHA), A, 
    as.double(LDA), B, as.double(LDB),
    as.double(BETA), C, as.double(LDC), as.logical(A.is.bm), 
    as.logical(B.is.bm), as.logical(C.is.bm), COFF)
}


# Vector addition and scaling
# Y := A * X  + Y
# A is a scalar double
# X is either a big.matrix or regular R matrix.
# Y is an optional matrix or vector of the same length as X.
# Returns a new matrix or big matrx with the same dimensions as X. If
# X is a dimension-less R vector, returns a column. Returned value type
# depends on the arguments and the value of the option
# options("bigalgebra.mixed_airthmetic_returns_R_matrix")[[1]].
daxpy = function(A=1, X, Y, LHS=1)
{
  mixed = FALSE
  X.is.bm = check_matrix(X,classes=c('big.matrix','matrix','vector','numeric'))
  # default to a column big matrix output
  M = length(X)
  L = M
  N = 1L
  D = dim(X)
  if(!is.null(D) && length(D)==2)
  {
    M = D[1]
    N = D[2]
  }
  Z = anon_matrix(M,N,val=0.0)
  if(!missing(Y))
  {
    # Check conformity of Y and duplicate
    if(length(Y)!=length(X)) stop("Lengths of X and Y must match")
    mixed = (X.is.bm != check_matrix(Y,classes=c('big.matrix','matrix','vector','numeric')))
    Z[] = Y[]
  }
  ans = daxpy_wrapper(as.double(L), as.double(A), X, Z, X.is.bm)
  if(mixed) return(ans[])
  ans
}

# Add a scalar to each element of a matrix
# Y := Y+SIGN*ALPHA 
dadd = function(Y, X, ALPHA)
{
  if (!is.numeric(ALPHA) || length(ALPHA) != 1)
    stop("ALPHA is not a scalar numeric value")

  Y.is.bm = check_matrix(Y)
  X.is.bm = check_matrix(X)
  
  D = dim(Y)
  M = D[1]
  N = D[2]
  
  Z <- anon_matrix(M,N,val=0.0)
  Z[] <- Y[]
  
  dadd_wrapper(as.double(ALPHA), Z, S, Y.is.bm)  
  return(S)
}


# Element-wise Matrix Multiplication
# C:= A * B
dgeemm = function(X, Y)
{
  if(class(X) != class(Y)){
    mixed = TRUE
  }else{
    mixed=FALSE
  }
  
  X.is.bm = check_matrix(X,classes=c('big.matrix','matrix','vector','numeric'))
  Y.is.bm = check_matrix(Y, classes=c('big.matrix', 'matrix', 'vector', 'numeric'))
    
  # Check conformity of matrices
  if(length(Y)!=length(X)) stop("Lengths of X and Y must match")
  
  # Dimensions
  D = dim(X)
  M = D[1]
  N = D[2]
  
  # create matrix for new values
  Z = anon_matrix(M,N,val=0.0)
  Z[] = Y[]

  ans = dgeemm_wrapper(X, Z, X.is.bm)
  
  if(mixed) return(ans[])
  return(ans)
}


# Element-wise Matrix Division
# C:= A / B
dgeemd = function(X, Y)
{
  if(class(X) != class(Y)){
    mixed = TRUE
  }else{
    mixed=FALSE
  }
  
  X.is.bm = is.big.matrix(X)
  
  # Check conformity of matrices
  if(length(Y)!=length(X)) stop("Lengths of X and Y must match")
  
  # Dimensions
  D = dim(X)
  M = D[1]
  N = D[2]
  
  # create matrix for new values
  Z = anon_matrix(M,N,val=0.0)
  Z[] = Y[]
  
  ans = dgeemd_wrapper(X, Z, X.is.bm)
  
  if(mixed) return(ans[])
  return(ans)
}

dgesmd = function(Y, ALPHA, ALPHA_LHS=0)
{
  if (!is.numeric(ALPHA) || length(ALPHA) != 1)
    stop("ALPHA is not a scalar numeric value")
  
  Y.is.bm = check_matrix(Y)
  
  # Dimensions
  D = dim(Y)
  M = D[1]
  N = D[2]
  
  # create matrix for new values
  Z = anon_matrix(M,N,val=0.0)
  Z[] = Y[]
  
  ans <- dgesmd_wrapper(as.double(ALPHA), Z, Y.is.bm, as.integer(ALPHA_LHS))
  return(ans)
}


# Matrix QR Decomposition
# Need to create a matrix that cannot be orthonormalized
# in order to properly test for exceptions when QR fails.
dgeqrf = function(A)
{
  A.is.bm = check_matrix(A,classes=c('big.matrix','matrix','vector','numeric'))
  
  # default to a column big matrix output
  D = dim(A)
  if(!is.null(D) && length(D)==2)
  {
    M = D[1]
    N = D[2]
  }
  
  # where A is m-by-n, produces an m-by-n upper triangular matrix R 
  # and an m-by-m unitary matrix Q so that A = Q*R.
  
  Q = anon_matrix(M,M,val=0.0)
  R <- anon_matrix(M,N,val=0.0)
  
#   ans = dgeqrf_wrapper(as.double(M), as.double(N), Y, LDA=as.double(M), 
#                        as.double(TAU), as.double(WORK), as.double(LWORK), as.double(INFO),
#                        as.logical(A.is.bm), as.logical(TAU.is.bm), as.logical(WORK.is.bm))

  ans = dgeqrf_wrapper(A, Q, R)
  return(ans)
}

# Cholesky factorization
# return 0 if successful, <0 if -i-th argument is invalid, > 0 if leading minor
# is not positive definite
# dpotrf=function(UPLO='U', N=NULL, A, LDA=NULL)
# {
#   if (is.null(N))
#   {
#     N = ncol(A)
#   }
#   if (is.null(LDA))
#   {
#     LDA = nrow(A)
#   }
#   A.is.bm = check_matrix(A)
#   INFO = 0
#   dpotrf_wrapper(as.character(UPLO), as.double(N), A, as.double(LDA),
#                  as.double(INFO), A.is.bm)
#   return(INFO)
# }

dpotrf=function(UPLO='U', N=NULL, A, LDA=NULL)
{
  if (is.null(N))
  {
    N = ncol(A)
  }
  if (is.null(LDA))
  {
    LDA = nrow(A)
  }
  
  A.is.bm = check_matrix(A)
  assert_is_positive_definite(A)
  
  
  R <- anon_matrix(N,N,val=0.0)
  R[] <- A[]
  
  INFO = 0
  dpotrf_wrapper(as.character(UPLO), as.double(N), R, as.double(LDA),
                 as.double(INFO), A.is.bm)
  return(R)
}


# transposition
transposeBM = function(X){
  
  assert_is_bigmatrix(X)
  
  # Dimensions
  D = dim(X)
  M = D[1]
  N = D[2]
  
  # create matrix for new values
  Y <- anon_matrix(N,M,val=0.0)

  ans <- t_wrapper(X,Y)
  
  return(ans)
}

# crosspord
bm_crossprod <- function(X, Y){
  if(nrow(X) != nrow(Y)){
    stop("matrices non-conformable")
  }
  
  Z <- anon_matrix(ncol(X), ncol(Y), val=0.0)
  cpp_bm_crossprod(X,Y,Z, 
                   is.big.matrix(X),
                   is.big.matrix(Y),
                   is.big.matrix(Z))
  
  return(Z)
}

# tcrosspord
bm_tcrossprod <- function(X, Y){
  if(ncol(X) != ncol(Y)){
    stop("matrices non-conformable")
  }
  
  Z <- anon_matrix(nrow(X), nrow(Y), val=0.0)
  cpp_bm_tcrossprod(X,Y,Z, 
                   is.big.matrix(X),
                   is.big.matrix(Y),
                   is.big.matrix(Z))
  
  return(Z)
}

# possible in-place function
# this will require a 'reshape' for the given
# big.matrix dimensions
# IPT_BM = function(X, low_mem=FALSE, inplace=FALSE){
#   
#   assert_is_bigmatrix(X)
#   assert_is_a_bool(low_mem)
#   assert_is_a_bool(inplace)
#   
#   transpose matrix
#     if(inplace){
#       ans <- t_inplace_wrapper(X,low_mem)
#     }else{
#       ans <- t_inplace_wrapper(Y,low_mem)
#     }
#   
#   return(ans)
# }


# Eigen
eigenBM <- function(X, only.values=TRUE){
  check_matrix(X, "big.matrix", "double")
  assert_is_a_bool(only.values)
  
  M = ncol(X)
  
  if(!only.values){
      EIG_VECS <- anon_matrix(M,M,val=0.0)  
  }else{
      EIG_VECS <- NULL
  }
  
  ret <- eigen_wrapper(X, EIG_VECS, only.values)
  return(ret)
}

# Power of matrix elements
# Y := POW(Y, B)
dgepow = function(Y, EXP)
{  
  if (!is.numeric(EXP) || length(EXP) != 1)
    stop("EXP is not a scalar numeric value")
    
  # Dimensions
  D = dim(Y)
  M = D[1]
  N = D[2]
  
  # create matrix for new values
  Z = anon_matrix(M,N,val=0.0)
  Z[] = Y[]
  
  ans <- dgepow_wrapper(EXP, Z)
  return(ans)
}


# Common log of matrix elements
# Y := LOG10(Y)
dgeclog = function(Y)
{
  # Dimensions
  D = dim(Y)
  M = D[1]
  N = D[2]
  
  # create matrix for new values
  Z = anon_matrix(M,N,val=0.0)
  Z[] = Y[]
  
  ans <- dgeclog_wrapper(Z)
  return(ans)
}


# Base log of matrix elements
# Y := LOG(Y, B)
dgelog = function(Y, BASE)
{
  # Dimensions
  D = dim(Y)
  M = D[1]
  N = D[2]
  
  # create matrix for new values
  Z = anon_matrix(M,N,val=0.0)
  Z[] = Y[]
  
  ans <- dgelog_wrapper(BASE, Z)
  return(ans)
}


# Exponential function of matrix elements
# Y := EXP(Y)
dgeexp = function(Y)
{
  # Dimensions
  D = dim(Y)
  M = D[1]
  N = D[2]
  
  # create matrix for new values
  Z = anon_matrix(M,N,val=0.0)
  Z[] = Y[]
  
  ans <- dgeexp_wrapper(Z)
  return(ans)
}


# Hyperbolic sine of matrix elements
# Y := SINH(Y)
dgesinh = function(Y)
{
  # Dimensions
  D = dim(Y)
  M = D[1]
  N = D[2]
  
  # create matrix for new values
  Z = anon_matrix(M,N,val=0.0)
  Z[] = Y[]
  
  ans <- dgesinh_wrapper(Z)
  return(ans)
}


# Hyperbolic cosine of matrix elements
# Y := COSH(Y)
dgecosh = function(Y)
{
  # Dimensions
  D = dim(Y)
  M = D[1]
  N = D[2]
  
  # create matrix for new values
  Z = anon_matrix(M,N,val=0.0)
  Z[] = Y[]
  
  ans <- dgecosh_wrapper(Z)
  return(ans)
}


# Hyperbolic tangent of matrix elements
# Y := TANH(Y)
dgetanh = function(Y)
{
  # Dimensions
  D = dim(Y)
  M = D[1]
  N = D[2]
  
  # create matrix for new values
  Z = anon_matrix(M,N,val=0.0)
  Z[] = Y[]
  
  ans <- dgetanh_wrapper(Z)
  return(ans)
}

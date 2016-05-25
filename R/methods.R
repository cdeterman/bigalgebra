#' @useDynLib bigalgebra
#' @import methods
#' @importFrom bigmemory typeof

# class unions when differentiating isn't important
setClassUnion("matrixORbigmatrix", c("matrix","big.matrix"))


#' @export
setMethod("%*%",signature(x="big.matrix", y="big.matrix"),
  function(x,y) 
  {
    dgemm(A=x, B=y)
  },
  valueClass="big.matrix"
)

#' @export
setMethod("%*%",signature(x="matrix", y="big.matrix"),
  function(x,y)
  {
    if(dim(x)[2] != dim(y)[1]) stop("non-conformant matrices")
    R = options("bigalgebra.mixed_airthmetic_returns_R_matrix")[[1]]
    if(!is.null(R) && R) return(dgemm(A=x, B=y, C=0))
    dgemm(A=x, B=y)
  },
  valueClass="matrix"
)

#' @export
setMethod("%*%",signature(x="big.matrix", y="matrix"),
  function(x,y) 
  {
    if(dim(x)[2] != dim(y)[1]) stop("non-conformant matrices")
    R = options("bigalgebra.mixed_airthmetic_returns_R_matrix")[[1]]
    if(!is.null(R) && R) return(dgemm(A=x, B=y, C=0))
    dgemm(A=x, B=y)
  },
  valueClass="matrix"
)

#' @export
setMethod("Arith",c(e1="big.matrix", e2="big.matrix"),
  function(e1,e2) 
  {
    op = .Generic[[1]]
    switch(op,
      `+` = daxpy(1.0,e1,e2),
      `-` = daxpy(-1.0,e2,e1),
      `*` = dgeemm(e1, e2),
      `/` = dgeemd(e1, e2),
      stop("Undefined operation")
    )
  }
)

#' @export
setMethod("Arith",c(e1="big.matrix", e2="matrix"),
  function(e1,e2) 
  {
    op = .Generic[[1]]
    switch(op,
      `+` = daxpy(1.0,e1,e2),
      `-` = daxpy(-1.0,e2,e1),
      `*` = dgeemm(e1, e2),
      `/` = dgeemd(e1, e2),
      stop("Undefined operation")
    )
  }
)

#' @export
setMethod("Arith",c(e1="matrix", e2="big.matrix"),
  function(e1,e2) 
  {
    op = .Generic[[1]]
    switch(op,
      `+` = daxpy(1.0,e1,e2),
      `-` = daxpy(-1.0,e2,e1),
      `*` = dgeemm(e2, e1),
      `/` = dgeemd(e1, e2),
      stop("Undefined operation")
    )
  }
)

#' @export
setMethod("Arith",c(e1="numeric", e2="big.matrix"),
  function(e1,e2) 
  {
    op = .Generic[[1]]
    if(length(e1)==1) {
      if (op=="*") 
        return(daxpy(e1,e2))
      return(switch(op,
        `+` = {
          e1 <- anon_matrix(nrow(e2), ncol(e2), val = e1)
          daxpy(1.0,e1,e2)
        },
        `-` = {
          e1 <- anon_matrix(nrow(e2), ncol(e2), val = e1)
          daxpy(-1.0,e2,e1)
        },
        `/` = dgesmd(e2, e1, 1),
        stop("Undefined operation")
      ))
    }
    stop("e1 is not a scalar")
  }
)

#' @export
setMethod("Arith",c(e1="big.matrix", e2="numeric"),
  function(e1,e2) 
  {
    op = .Generic[[1]]
    if(length(e2)==1) {
      if( op=="*") 
        return(daxpy(e2,e1))
      return(switch(op,
        `+` = {
          e2 <- anon_matrix(nrow(e1), ncol(e1), val = e2)
          daxpy(1.0,e1,e2)
        },
        `-` = {
          e2 <- anon_matrix(nrow(e1), ncol(e1), val = e2)
          daxpy(-1.0,e2,e1)
        },
        `^` = dgepow(e1, e2),
        `/` = dgesmd(e1, e2, 0),
        stop("Undefined operation")
      ))
    }
    stop("e2 is not a scalar")
  }
)


#' @export
setMethod("Arith", c(e1="big.matrix", e2="missing"),
          function(e1, e2)
          {
            op = .Generic[[1]]
            switch(op,
                   `-` = unary_axpy(e1),
                   stop("undefined operation")
            )
          },
          valueClass = "gpuMatrix"
)

#' @export
setMethod("Math", c(x="big.matrix"),
          function(x)
          {
            op = .Generic[[1]]
            switch(op,
                   `sqrt` = dgepow(x, 0.5),
                   `exp` = dgeexp(x),
                   `log10` = dgeclog(x),
                   `sinh` = dgesinh(x),
                   `cosh` = dgecosh(x),
                   `tanh` = dgetanh(x),
                   stop("Undefined operation")
            )
          }
)

#' @export
setMethod("log", signature(x="big.matrix"),
          function(x, base=exp(1))
          {
            dgelog(x, base)
          }
)

#' @export
setMethod("qr", c(x="big.matrix"),
          function(x, ...)
            {
              dgeqrf(x)
            })

#' @export
setMethod("chol", c(x="big.matrix"),
          function(x)
          {
            dpotrf(A=x)
          })

#' @export
setMethod("eigen", c(x="big.matrix"),
          function(x, only.values=TRUE)
            {
              eigenBM(x, only.values)
            })

#' @export
t.big.matrix <- 
  function(x){
    transposeBM(x)
  }

#' @title big.matrix crossproduct
#' @description Return the matrix cross-product of two conformable
#' big.matrix objects.  This is equivalent to t(x) %*% y (crossprod)
#' or x %*% t(t) (tcrossprod) but faster as no data transfer between
#' device and host is required.
#' @param x A big.matrix
#' @param y A big.matrix
#' @return A big.matrix
#' @author Charles Determan Jr.
#' @docType methods
#' @rdname big.matrix-crossprod
#' @aliases crossprod,big.matrix
#' @export
setMethod("crossprod",
          signature(x = "big.matrix", y = "missing"),
          function(x, y){
            bm_crossprod(x, x)
          })

#' @rdname big.matrix-crossprod
#' @export
setMethod("crossprod",
          signature(x = "big.matrix", y = "big.matrix"),
          function(x, y){
            bm_crossprod(x, y)
          })

#' @rdname big.matrix-crossprod
setMethod("tcrossprod",
          signature(x = "big.matrix", y = "missing"),
          function(x, y){
            bm_tcrossprod(x, x)
          })


#' @rdname big.matrix-crossprod
#' @export
setMethod("tcrossprod",
          signature(x = "big.matrix", y = "big.matrix"),
          function(x, y){
            bm_tcrossprod(x, y)
          })

#' @title all.equal
#' @description all.equal method for comparing elements between big.matrix and
#' matrix objects.
#' @param target A \code{matrix} or \code{big.matrix} object
#' @param current A \code{matrix} or \code{big.matrix} object
#' @param tolerance tolerance for rounding error
#' @return Either \code{TRUE} or \code{FALSE} indicating if all elements are
#' equivlanet.
#' @export
setMethod("all.equal", c(target = "matrixORbigmatrix", current="matrixORbigmatrix"),
          function(target,current,tolerance=.Machine$double.eps^0.5){
            if(!all(dim(target) == dim(current))){
              return(FALSE)
            }
            X_isBM <- check_matrix(target)
            Y_isBM <- check_matrix(current)
            
            all_equal_cpp(target, current, X_isBM, Y_isBM, tolerance)
          }
)


#' @export
setMethod("isSquare", c(object="big.matrix"),
          function(object){
            if(!is.big.matrix(object)){
              stop("argument is not a big.matrix")
            }
            return(nrow(object) == ncol(object))
          }
)

#' @export
setMethod("isSquare", c(object="matrix"),
          function(object){
            if(!is.matrix(object)){
              stop("argument is not a matrix")
            }
            return(nrow(object) == ncol(object))
          }
)

#' @export
setMethod("isDiagonal", c(object="big.matrix"),
          function(object){
              if(!isSquare(object)){
                return(FALSE)
              }
              
              Y.is.bm = check_matrix(object)
              cpp_isDiagonal(object, Y.is.bm)
          })

#' @export
setMethod("isDiagonal", c(object="matrix"),
          function(object){
            if(!isSquare(object)){
              return(FALSE)
            }
            
            Y.is.bm = check_matrix(object)
            cpp_isDiagonal(object, Y.is.bm)
          })

#' @export
setMethod("isTriangular", c(object="big.matrix"),
          function(object, upper=NA){
            if(!isSquare(object)){
              return(FALSE)
            }
            
            Y.is.bm = check_matrix(object)
            cpp_isTriangular(object, Y.is.bm, upper)
          })

#' @export
setMethod("isTriangular", c(object="matrix"),
          function(object, upper=NA){
            if(!isSquare(object)){
              return(FALSE)
            }
            
            Y.is.bm = check_matrix(object)
            cpp_isTriangular(object, Y.is.bm, upper)
          })

#' @export
setMethod("isSymmetric", signature(object="big.matrix"),
          function(object, tol = 100 * .Machine$double.eps, ...){
            
            if (!is.big.matrix(object)) {
              stop("argument x is not a big.matrix")
            }
            
            check_matrix(object, classes="big.matrix", types=c("integer", "double"))
            
            if (!isSquare(object)) 
              stop("argument object is not a square numeric matrix")
            
            return(all.equal(object, t(object)))
          }
)

#' @export
setMethod("isPositiveDefinite", signature(object="big.matrix"),
          function(object, tol = 100 * .Machine$double.eps){
            
            if (!is.big.matrix(object)) {
              stop("argument x is not a big.matrix")
            }
            
            check_matrix(object, classes="big.matrix", types=c("integer", "double"))
            
            if (!isSquare(object)) 
                stop("argument object is not a square numeric matrix")
            if (!isSymmetric(object, tol)) 
                stop("argument x is not a symmetric matrix")
            
            eigenvalues <- eigen(object, only.values = TRUE)$values
            
            n <- nrow(object)
            for (i in 1:n) {
              if (abs(eigenvalues[i]) < tol) {
                eigenvalues[i] <- 0
              }
            }
            
            if (any(eigenvalues <= 0)) {
              return(FALSE)
            }
            return(TRUE)
          }
)

#' @export
setMethod("isPositiveDefinite", signature(object="matrix"),
          function(object, tol = 100 * .Machine$double.eps){
            
            if (!is.matrix(object)) {
              stop("argument x is not a matrix")
            }
            
            check_matrix(object, classes="matrix", types=c("integer", "double"))
            
            if (!isSquare(object)) 
              stop("argument object is not a square numeric matrix")
            if (!isSymmetric(object, tol)) 
              stop("argument x is not a symmetric matrix")
            
            eigenvalues <- eigen(object, only.values = TRUE)$values
            
            n <- nrow(object)
            for (i in 1:n) {
              if (abs(eigenvalues[i]) < tol) {
                eigenvalues[i] <- 0
              }
            }
            
            if (any(eigenvalues <= 0)) {
              return(FALSE)
            }
            return(TRUE)
          }
)

#' @export
setMethod("isPositiveSemiDefinite", signature(object="big.matrix"),
          function(object, tol = 100 * .Machine$double.eps){
            
            if (!is.big.matrix(object)) {
              stop("argument x is not a big.matrix")
            }
            
            check_matrix(object, classes="big.matrix", types=c("integer", "double"))
            
            if (!isSquare(object)) 
              stop("argument object is not a square numeric matrix")
            if (!isSymmetric(object, tol)) 
              stop("argument x is not a symmetric matrix")
            
            eigenvalues <- eigen(object, only.values = TRUE)$values
            
            n <- nrow(object)
            for (i in 1:n) {
              if (abs(eigenvalues[i]) < tol) {
                eigenvalues[i] <- 0
              }
            }
            
            if (any(eigenvalues < 0)) {
              return(FALSE)
            }
            return(TRUE)
          }
)

#' @export
setMethod("isPositiveSemiDefinite", signature(object="matrix"),
          function(object, tol = 100 * .Machine$double.eps){
            
            if (!is.matrix(object)) {
              stop("argument x is not a matrix")
            }
            
            check_matrix(object, classes="matrix", types=c("integer", "double"))
            
            if (!isSquare(object)) 
              stop("argument object is not a square numeric matrix")
            if (!isSymmetric(object, tol)) 
              stop("argument x is not a symmetric matrix")
            
            eigenvalues <- eigen(object, only.values = TRUE)$values
            
            n <- nrow(object)
            for (i in 1:n) {
              if (abs(eigenvalues[i]) < tol) {
                eigenvalues[i] <- 0
              }
            }
            
            if (any(eigenvalues < 0)) {
              return(FALSE)
            }
            return(TRUE)
          }
)

#' @export
setMethod("isNegativeDefinite", signature(object="big.matrix"),
          function(object, tol = 100 * .Machine$double.eps){
            
            if (!is.big.matrix(object)) {
              stop("argument x is not a big.matrix")
            }
            
            check_matrix(object, classes="big.matrix", types=c("integer", "double"))
            
            if (!isSquare(object)) 
              stop("argument object is not a square numeric matrix")
            if (!isSymmetric(object, tol)) 
              stop("argument x is not a symmetric matrix")
            
            eigenvalues <- eigen(object, only.values = TRUE)$values
            
            n <- nrow(object)
            for (i in 1:n) {
              if (abs(eigenvalues[i]) < tol) {
                eigenvalues[i] <- 0
              }
            }
            
            if (any(eigenvalues >= 0)) {
              return(FALSE)
            }
            return(TRUE)
          }
)

#' @export
setMethod("isNegativeDefinite", signature(object="matrix"),
          function(object, tol = 100 * .Machine$double.eps){
            
            if (!is.matrix(object)) {
              stop("argument x is not a matrix")
            }
            
            check_matrix(object, classes="matrix", types=c("integer", "double"))
            
            if (!isSquare(object)) 
              stop("argument object is not a square numeric matrix")
            if (!isSymmetric(object, tol)) 
              stop("argument x is not a symmetric matrix")
            
            eigenvalues <- eigen(object, only.values = TRUE)$values
            
            n <- nrow(object)
            for (i in 1:n) {
              if (abs(eigenvalues[i]) < tol) {
                eigenvalues[i] <- 0
              }
            }
            
            if (any(eigenvalues >= 0)) {
              return(FALSE)
            }
            return(TRUE)
          }
)

#' @export
setMethod("isNegativeSemiDefinite", signature(object="big.matrix"),
          function(object, tol = 100 * .Machine$double.eps){
            
            if (!is.big.matrix(object)) {
              stop("argument x is not a big.matrix")
            }
            
            check_matrix(object, classes="big.matrix", types=c("integer", "double"))
            
            if (!isSquare(object)) 
              stop("argument object is not a square numeric matrix")
            if (!isSymmetric(object, tol)) 
              stop("argument x is not a symmetric matrix")
            
            eigenvalues <- eigen(object, only.values = TRUE)$values
            
            n <- nrow(object)
            for (i in 1:n) {
              if (abs(eigenvalues[i]) < tol) {
                eigenvalues[i] <- 0
              }
            }
            
            if (any(eigenvalues > 0)) {
              return(FALSE)
            }
            return(TRUE)
          }
)

#' @export
setMethod("isNegativeSemiDefinite", signature(object="matrix"),
          function(object, tol = 100 * .Machine$double.eps){
            
            if (!is.matrix(object)) {
              stop("argument x is not a matrix")
            }
            
            check_matrix(object, classes="matrix", types=c("integer", "double"))
            
            if (!isSquare(object)) 
              stop("argument object is not a square numeric matrix")
            if (!isSymmetric(object, tol)) 
              stop("argument x is not a symmetric matrix")
            
            eigenvalues <- eigen(object, only.values = TRUE)$values
            
            n <- nrow(object)
            for (i in 1:n) {
              if (abs(eigenvalues[i]) < tol) {
                eigenvalues[i] <- 0
              }
            }
            
            if (any(eigenvalues > 0)) {
              return(FALSE)
            }
            return(TRUE)
          }
)

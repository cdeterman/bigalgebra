### additional classes of big.matrix objects
# These are a work in progress at the moment with no current
# implementation.  Similar to the Matrix package, these classes
# will be useful for providing more efficient analysis of special
# matrix types.


setClass("symmetricBigMatrix", contains="big.matrix",
         validity = function(object) isSymmetric(object))
setClass("triangularBigMatrix", contains="big.matrix",
         validity = function(object) isTriangular(object))

### compact classes, similar to Matrix package
# double, symmetric, big.matrix
setClass('dsyBigMatrix', contains='symmetricBigMatrix')
# double, triangular, big.matrix
setClass('dtriBigMatrix', contains='triangularBigMatrix')

# # integer, symmetric, big.matrix
# setClass('isyBigMatrix', contains='symmetricBigMatrix')
# # integer, triangular, big.matrix
# setClass('itriBigMatrix', contains='triangularBigMatrix')
# 
# # logical, symmetric, big.matrix
# setClass('lsyBigMatrix', contains='symmetricBigMatrix')
# # logical, triangular, big.matrix
# setClass('ltriBigMatrix', contains='triangularBigMatrix')
# 
# # complex, symmetric, big.matrix
# setClass('zsyBigMatrix', contains='symmetricBigMatrix')
# # complex, triangular, big.matrix
# setClass('ztriBigMatrix', contains='triangularBigMatrix')

### alternate verbose style???
# setClass('big.symmetric.matrix', contains='big.matrix')
# setClass('big.triangular.matrix', contains='big.matrix')
# setClass('big.diagonal.matrix', contains='big.matrix')
# setClass('big.positive.definite.matrix', contains='big.matrix')

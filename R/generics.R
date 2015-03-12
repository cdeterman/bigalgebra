# Generics for various functions

# These matrix type functions seem like they should be 
# in bigmemory but internals depend on bigalgebra so they 
# will remain here for now
setGeneric("isSquare", function(object){
  standardGeneric("isSquare")
})

setGeneric("isDiagonal", function(object){
  standardGeneric("isDiagonal")
})

setGeneric("isTriangular", function(object, ...){
  standardGeneric("isTriangular")
})

setGeneric("isPositiveDefinite", function(object, ...){
  standardGeneric("isPositiveDefinite")
})

setGeneric("isPositiveSemiDefinite", function(object, ...){
  standardGeneric("isPositiveSemiDefinite")
})

setGeneric("isNegativeDefinite", function(object, ...){
  standardGeneric("isNegativeDefinite")
})

setGeneric("isNegativeSemiDefinite", function(object, ...){
  standardGeneric("isNegativeSemiDefinite")
})

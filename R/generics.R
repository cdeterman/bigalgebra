# Generics for various functions

# These matrix type functions seem like they should be 
# in bigmemory but internals depend on bigalgebra so they 
# will remain here for now
setGeneric("isSquare", function(object){
  standardGeneric("isSquare")
})

setGeneric("isPositiveDefinite", function(object, ...){
  standardGeneric("isPositiveDefinite")
})

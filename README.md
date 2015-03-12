bigalgebra
==========

### WISH List
1. Use of additional BLAS libraries including ATLAS and OpenBLAS
  a. This is likely solved by utilizing RcppArmadillo which automatically
  will interface with whatever BLAS implementation R is pointing towards.
2. Option to interface with GPU (CUDA, OpenCL, ArrayFire?), likely will need 
to be a separate package to avoid compiling issues.

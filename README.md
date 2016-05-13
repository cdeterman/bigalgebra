
bigmemory
=========

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/bigalgebra)](http://cran.r-project.org/package=bigalgebra)
[![Build Status](https://travis-ci.org/kaneplusplus/bigalgebra.png)](https://travis-ci.org/kaneplusplus/bigalgebra)
[![Build status](https://ci.appveyor.com/api/projects/status/uybv01gtro3xl58y/branch/master?svg=true)](https://ci.appveyor.com/project/kaneplusplus/bigalgebra/branch/master)
[![Coverage Status](https://coveralls.io/repos/kaneplusplus/bigalgebra/badge.svg?branch=master&service=github)](https://coveralls.io/github/kaneplusplus/bigalgebra?branch=master)


### WISH List
1. Use of additional BLAS libraries including ATLAS and OpenBLAS
  a. This is likely solved by utilizing RcppArmadillo which automatically
  will interface with whatever BLAS implementation R is pointing towards.
2. Option to interface with GPU (CUDA, OpenCL, ArrayFire?), likely will need 
to be a separate package to avoid compiling issues.
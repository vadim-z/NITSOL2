Fork of NITSOL project by H.F.Walker:
Taken from [here](http://users.wpi.edu/~walker/NITSOL/)

Changes:
* classical modified Gram-Schmidt GMRES is implemented to complement simpler GMRES
* finite-difference precision can be specified to handle RHS which is
evaluated with lower accuracy
* common blocks are eliminated, all inputs and outputs are passed as parameters
* example to test Krylov solvers is added
* pieces of BLAS and LAPACK removed, user should link to external BLAS and LAPACK
* example CMake build file is provided

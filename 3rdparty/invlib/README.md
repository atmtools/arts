# invlib

invlib is a C++ template library that provides a generic implementation of maximum
a-posteriori (MAP) estimators for inverse problems. The aim of the library is to
provide a formulation of MAP estimators independent of the data types used for the
underlying matrix algebra. This is required to be able to efficiently handle the
many different problem sizes that occur in the retrieval of atmospheric remote
sensing data.

invlib is currently under development and will be integrated into the open-source
atmospheric radiative transfer code [ARTS](http://www.radiativetransfer.org/).

## Current Features:

- **Generic Formulation**: Supports arbitrary matrix and vector types, as well
  as separate types for Jacobian and covariance matrices.
- **Standard**, **n-form** and **m-form** formulation of the MAP estimators.
- **Optimization Methods**: Gauss-Newton, Levenberg-Marquardt
- **Indirect solver**: Provides a conjugate gradient solver for the subproblem
    of the Gauss-Newton or Levenberg-Marquardt subproblem.
- **Optimized symbolic matrix algebra**: The symbolic algebra delays computation
    of matrix arithmetic operations until it gets converted to the result type.
    This avoids the computation of expensive matrix-matrix products if the result
    of the expression is actually a vector.

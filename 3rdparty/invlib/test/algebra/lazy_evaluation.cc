#include "invlib/interfaces/arts_wrapper.h"
#include "invlib/algebra.h"
#include "invlib/algebra/solvers.h"
#include "invlib/algebra/precision_matrix.h"
#include "utility.h"

#include <stdio.h>
#include <iostream>

template <typename T> void foo(T a);

int main( int argc, const char** argv )
{

    using Vector   = invlib::Vector<ArtsVector>;
    using Matrix   = invlib::Matrix<ArtsMatrix<Matrix>, Vector>;
    using Identity = invlib::MatrixIdentity<double, Matrix>;
    using PrecisionMatrix = invlib::PrecisionMatrix<Matrix>;
    const unsigned int n = 5;

    Matrix A, B, C, D, E, F, G, H;
    Identity I;
    Vector v, w;

    A.resize(n,n); B.resize(n,n); C.resize(n,n); D.resize(n,n);
    E.resize(n,n); F.resize(n,n); G.resize(n,n); H.resize(n,n);

    v.resize(n); w.resize(n);

    A.resize(1,5);
    B.resize(5,5);

    A(0,0) = 1.0;
    A(0,1) = 1.0;
    A(0,2) = 1.0;

    B(0,0) = 1.0;
    B(1,1) = 1.0;
    B(2,2) = 1.0;
    B(0,1) = 1.0;
    B(0,2) = 1.0;

    C = transp(B) * transp(A);

    std::cout << A << std::endl << " --- " << std::endl;
    std::cout << B << std::endl << " --- " << std::endl;
    std::cout << C << std::endl << " --- " << std::endl;

    ::Matrix X(1,5), Y(5,5), Z(5,1);
    X(0,0) = 1.0;
    X(0,1) = 1.0;
    X(0,2) = 1.0;
    ::Matrix XX(transpose(X));

    Y(0,0) = 1.0;
    Y(1,1) = 1.0;
    Y(2,2) = 1.0;
    Y(0,1) = 1.0;
    Y(0,2) = 1.0;

    mult(Z, transpose(Y), XX);

    std::cout << X << std::endl << " --- " << std::endl;
    std::cout << Y << std::endl << " --- " << std::endl;
    std::cout << Z << std::endl << " --- " << std::endl;
    //printf("H(0,0) = %f \n v(0) = %f \n", maximum_error(w,v), w(0));
}

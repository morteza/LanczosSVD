#ifndef EIGENVALUEDECOMPOSITION_H
#define EIGENVALUEDECOMPOSITION_H

#include "Matrix.h"
#include "Vector.h"

class EigenValueDecomposition
{
public:
    EigenValueDecomposition(const Matrix &input);

    //! input matrix dimension
    long n;

    //! storage of eigenvalues (Real parts)
    double *d;
    //! storage of eigenvalues (Imaginary parts)
    double *e;

    //! storage of eigenvectors
    double **V;

private:
    //! Symmetric Householder reduction to tridiagonal form
    void triDiagonalize();

    //! Symmetric tridiagonal QL decomposition algorithm
    void QLDecomposition();
    //! private method for sqrt(a^2 + b^2) without under/overflow
    double hypot(double a, double b);
};

#endif // EIGENVALUEDECOMPOSITION_H

#ifndef LINEARSYSTEM_H
#define LINEARSYSTEM_H

#include "activeRodMesh.h"
#include "fvCFD.H"
#include <limits>

namespace Foam
{
    void linearSolve
    (
        const gismo::gsMatrix<scalar>& A,
        gismo::gsMatrix<scalar>& x,
        const gismo::gsMatrix<scalar>& b,
        scalar tol=1e-10
    );
    
    void linearSolve
    (
        const gismo::gsMatrix<scalar>& A,
        gismo::gsMatrix<scalar>& x,
        const gismo::gsMatrix<scalar>& b,
        const gismo::gsMatrix<scalar>& x_initial,
        scalar tol=1e-10
    );
    
    void linearSolve_ConjugateGradient
    (
        const gismo::gsMatrix<scalar>& A,
        gismo::gsMatrix<scalar>& x,
        const gismo::gsMatrix<scalar>& b,
        const gismo::gsMatrix<scalar>& x_initial,
        scalar tol=1e-10
    );
    
    void linearSolve_GaussSeidel
    (
        const gismo::gsMatrix<scalar>& A,
        gismo::gsMatrix<scalar>& x,
        const gismo::gsMatrix<scalar>& b,
        const gismo::gsMatrix<scalar>& x_initial,
        scalar tol=1e-10
    );
}

#endif

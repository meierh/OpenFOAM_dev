#ifndef LINEARSYSTEM_H
#define LINEARSYSTEM_H

#include "activeRodMesh.h"
#include "fvCFD.H"

namespace Foam
{
    void linearSolve
    (
        gismo::gsMatrix<scalar>& A,
        gismo::gsMatrix<scalar>& x,
        gismo::gsMatrix<scalar>& b,
        scalar tol=1e-10
    );
}

#endif

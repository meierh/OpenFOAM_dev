#include "LinearSystem.H"

void Foam::linearSolve
(
    gismo::gsMatrix<scalar>& A,
    gismo::gsMatrix<scalar>& x,
    gismo::gsMatrix<scalar>& b,
    scalar tol
)
{
    if(A.cols()!=A.rows())
        FatalErrorInFunction<<"Matrix must be square"<<exit(FatalError);
    if(A.cols()!=b.rows())
        FatalErrorInFunction<<"Matrix and rhs dimensions dont match"<<exit(FatalError);
    if(b.cols()!=1)
        FatalErrorInFunction<<"rhs must be a vector"<<exit(FatalError);
    
    x = gismo::gsMatrix<scalar>(b.rows(),1);
    
    
    
}

#include "LinearSystem.H"

void Foam::linearSolve
(
    const gismo::gsMatrix<scalar>& A,
    gismo::gsMatrix<scalar>& x,
    const gismo::gsMatrix<scalar>& b,
    scalar tol
)
{    
    gismo::gsMatrix<scalar> x_initial(b.rows(),1);
    for(label i=0; i<b.rows(); i++)
        x_initial(i,0) = 0;
    linearSolve(A,x,b,x_initial,tol);
}

void Foam::linearSolve
(
    const gismo::gsMatrix<scalar>& A,
    gismo::gsMatrix<scalar>& x,
    const gismo::gsMatrix<scalar>& b,
    const gismo::gsMatrix<scalar>& x_initial,
    scalar tol
)
{
    if(A.cols()!=A.rows())
        FatalErrorInFunction<<"Matrix must be square"<<exit(FatalError);
    if(A.cols()!=b.rows())
        FatalErrorInFunction<<"Matrix and rhs dimensions dont match"<<exit(FatalError);
    if(b.cols()!=1)
        FatalErrorInFunction<<"rhs must be a vector"<<exit(FatalError);
    
    x = x_initial;
        
    gismo::gsConjugateGradient Solver(A);
    Solver.setTolerance(tol);
    Solver.solve(b,x);
    
    gismo::gsMatrix<scalar> resid = A*x-b;
    scalar residNormL2=0;
    for(label i=0;i<resid.rows();i++)
        residNormL2+=(resid(i,0)*resid(i,0));
    residNormL2 = std::sqrt(residNormL2);
    
    scalar xNormL2=0;
    for(label i=0;i<x.rows();i++)
        xNormL2+=(x(i,0)*x(i,0));
    xNormL2 = std::sqrt(xNormL2);
    
    /*
    std::cout<<std::endl;
    std::cout<<"residNormL2:"<<residNormL2<<std::endl;
    std::cout<<"xNormL2:"<<xNormL2<<std::endl;
    std::cout<<"residNormL2/xNormL2:"<<residNormL2/xNormL2<<std::endl;
    std::cout<<"Solver.error():"<<Solver.error()<<std::endl;
    std::cout<<std::endl;
    */

    if(Solver.error()>1e-8 || residNormL2>1e-8)
    {
        std::cout<<"A:"<<std::endl<<A<<std::endl;
        std::cout<<"x:"<<std::endl<<x<<std::endl;
        std::cout<<"b:"<<std::endl<<b<<std::endl;
        std::cout<<"A*x:"<<std::endl<<A*x<<std::endl;
        std::cout<<"residNormL2:"<<residNormL2<<std::endl;
        std::cout<<"Solver.iterations():"<<Solver.iterations()<<std::endl;
        std::cout<<"Solver.error():"<<Solver.error()<<std::endl;
        
        linearSolve_gaussSeidel(A,x,b,x);
        
        std::cout<<"A:"<<std::endl<<A<<std::endl;
        std::cout<<"x:"<<std::endl<<x<<std::endl;
        std::cout<<"b:"<<std::endl<<b<<std::endl;
        std::cout<<"A*x:"<<std::endl<<A*x<<std::endl;
        
        FatalErrorInFunction<<"Temp Stop"<<exit(FatalError);

        FatalErrorInFunction<<"Failed computing weights"<<exit(FatalError);
    }
}
    
void Foam::linearSolve_ConjugateGradient
(
    const gismo::gsMatrix<scalar>& A,
    gismo::gsMatrix<scalar>& x,
    const gismo::gsMatrix<scalar>& b,
    const gismo::gsMatrix<scalar>& x_initial,
    scalar tol
)
{
    
}
    
void Foam::linearSolve_gaussSeidel
(
    const gismo::gsMatrix<scalar>& A,
    gismo::gsMatrix<scalar>& x,
    const gismo::gsMatrix<scalar>& b,
    const gismo::gsMatrix<scalar>& x_initial,
    scalar tol
)
{
    scalar error = std::numeric_limits<scalar>::max();
    x = x_initial;
    label iterationCounter=0;
    
    while(error>tol)
    {
        for(label k=0; k<A.rows(); k++)
        {
            scalar sumk_Ax = 0;
            for(label i=0; i<A.cols(); i++)
            {
                sumk_Ax += A(k,i)*x(i,0);
            }
            x(k) = (b(k,0)-sumk_Ax)/A(k,k);
        }
        
        gismo::gsMatrix<scalar> resid = A*x-b;
        scalar residNormL2=0;
        for(label i=0;i<resid.rows();i++)
            residNormL2+=(resid(i,0)*resid(i,0));
        residNormL2 = std::sqrt(residNormL2);
        
        scalar xNormL2=0;
        for(label i=0;i<x.rows();i++)
            xNormL2+=(x(i,0)*x(i,0));
        xNormL2 = std::sqrt(xNormL2);
        
        error = residNormL2/xNormL2;
        
        iterationCounter++;
        if(iterationCounter>=100000)
            Pout<<"Warning: GaussSeidel Iteration high"<<Foam::endl;
    }
}

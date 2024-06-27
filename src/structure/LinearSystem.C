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
        
        linearSolve_GaussSeidel(A,x,b,x);
        
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
    
void Foam::linearSolve_GaussSeidel
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
            for(label i=0; i<k; i++)
            {
                sumk_Ax += A(k,i)*x(i,0);
            }
            for(label i=k+1; i<A.cols(); i++)
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
            Pout<<"Warning: GaussSeidel Iteration high | error:"<<error<<Foam::endl;
    }
}

void Foam::CSR_Matrix_par::addRow(List<scalar> row)
{
    if(row.size()!=globalCols)
        FatalErrorInFunction<<"Invalid row size"<<exit(FatalError);
    
    DynamicList<std::pair<scalar,label>> rowData;
    for(label col=0; col<row.size(); col++)
    {
        if(row[col]!=0)
        {
            rowData.append({row[col],col});
        }
    }
    addRow(rowData);
}

void Foam::CSR_Matrix_par::addRow(List<std::pair<scalar,label>> row)
{
    if(writtenRow<localRows)
        writtenRow++;
    else
        FatalErrorInFunction<<"Out of range row"<<exit(FatalError);
    
    Row_Index.append(lastWrittenRowIndex);
    lastWrittenRowIndex+=row.size();
    for(label i=0; i<row.size(); i++)
    {
        std::pair<scalar,label>& entry = row[i];
        scalar value = entry.first;
        label col = entry.second;
        if(col>=globalCols)
            FatalErrorInFunction<<"Col out of range"<<exit(FatalError);
        V.append(value);
        Col_Index.append(col);
    }
}

Vector_par Foam::CSR_Matrix_par::operator*
(
    const Vector_par& vec
) const
{
    checkCompatible(vec);
    
    List<List<scalar>> globalVecs(Pstream::nProcs());    
    std::pair<label,label> vecInfo = vec.getLocalSize();
    label localRowStart = vecInfo.first;
    label localRows = vecInfo.second;
    globalVecs[Pstream::myProcNo()] = List<scalar>(localRows);
    for(label i=0; i<localRows; i++)
    {
        globalVecs[Pstream::myProcNo()][i] = vec[i];
    }
    Pstream::gatherList(globalVecs);
    Pstream::scatterList(globalVecs);
    
    DynamicList<scalar> globalVec;
    for(label proc=0; proc<Pstream::myProcNo()-1; proc++)
    {
        const List<scalar>& procVec = globalVecs[proc];    
        for(label ind=0; ind<procVec.size(); ind++)
        {
            globalVec.append(procVec[ind]);
        }
    }
    if(globalVec.size()!=vec.getLocalSize().first)
        FatalErrorInFunction<<"Global to local size mismatch"<<exit(FatalError);
    for(label proc=Pstream::myProcNo(); proc<Pstream::nProcs(); proc++)
    {
        const List<scalar>& procVec = globalVecs[proc];    
        for(label ind=0; ind<procVec.size(); ind++)
        {
            globalVec.append(procVec[ind]);
        }
    }
    if(globalVec.size()!=vec.getGlobalSize())
        FatalErrorInFunction<<"Global size mismatch"<<exit(FatalError);
    
    Vector_par result(this->localRows,this->localRowStart,this->globalRows);
    for(label k=0; k<this->localRows; k++)
    {
        label rowStartInd = Row_Index[k];
        label rowEndInd;
        if(k<this->localRows-1)
            rowEndInd = Row_Index[k+1];
        else
            rowEndInd = V.size();
        
        scalar value = 0;
        for(label ind=rowStartInd; ind<rowEndInd; ind++)
        {
            label col = Col_Index[ind];
            value += V[ind]*globalVec[col];
        }
        result[k] = value;
    }
    return result;
}

CSR_Matrix_par Foam::CSR_Matrix_par::operator* 
(
    const CSR_Matrix_par& mat
) const
{
    FatalErrorInFunction<<"Not yet implemented"<<exit(FatalError);
}

scalar Foam::CSR_Matrix_par::getDiagonal
(
    label localRow
) const
{
    label rowStartInd = Row_Index[localRow];
    label rowEndInd;
    if(localRow<this->localRows-1)
        rowEndInd = Row_Index[localRow+1];
    else
        rowEndInd = V.size();
    for(label ind=rowStartInd; ind<rowEndInd; ind++)
    {
        label col = Col_Index[ind];
        if(col==localRowStart+localRow)
            return V[ind];
    }
    return 0;    
}

scalar Foam::CSR_Matrix_par::evalOnDiagonal
(
    const Vector_par& vec,
    label localRow
) const
{
    checkCompatible(vec);
    
    scalar diag_A = getDiagonal(localRow);
    if(diag_A==0)
        return 0;
    else
        return diag_A*vec[localRow];
}

void Foam::CSR_Matrix_par::checkCompatible
(
    const Vector_par& vec
) const
{
    if(globalCols!=vec.getGlobalSize())
        FatalErrorInFunction<<"Matrix Vector mismatch"<<exit(FatalError);
}

void Foam::CSR_Matrix_par::checkCompatible
(
    const CSR_Matrix_par& mat
) const
{
    if(globalCols!=mat.globalRows)
        FatalErrorInFunction<<"Matrix matrix mismatch"<<exit(FatalError);
}

scalar& Foam::Vector_par::operator[]
(
    int idx
)
{
    return V[idx];
}

scalar Foam::Vector_par::operator[]
(
    int idx
) const
{
    return V[idx];
}

scalar Foam::Vector_par::operator&
(
    const Vector_par& vec
) const
{
    checkCompatible(vec);
    
    scalar localValue = 0;
    for(label i=0; i<localRows; i++)
        localValue += (*this)[i]*vec[i];
    
    struct scalar_plus
    {
        scalar operator()(const scalar a, const scalar b) const
        {
            return a+b;
        }
    };
    scalar_plus op;
    Pstream::gather(localValue,op);
    Pstream::scatter(localValue);
    return localValue;
}

Vector_par Foam::Vector_par::operator*
(
    scalar scale
) const
{
    Vector_par result(*this);
    for(label i=0; i<localRows; i++)
        result[i] = (*this)[i]*scale;
    return result;
}

Vector_par Foam::Vector_par::operator+
(
    const Vector_par& vec
) const
{
    checkCompatible(vec);
    Vector_par result(vec);
    for(label i=0; i<localRows; i++)
        result[i] = (*this)[i]+vec[i];
    return result;
}

Vector_par Foam::Vector_par::operator-
(
    const Vector_par& vec
) const
{
    checkCompatible(vec);
    Vector_par result(vec);
    for(label i=0; i<localRows; i++)
        result[i] = (*this)[i]-vec[i];
    return result;
}

scalar Foam::Vector_par::norm2() const
{
    return std::sqrt((*this)&(*this));
}

void Foam::Vector_par::scale
(
    scalar scale
)
{
    for(label i=0; i<localRows; i++)
        (*this)[i] *= scale;
}

void Foam::Vector_par::fill
(
    scalar value
)
{
    for(label i=0; i<localRows; i++)
        (*this)[i] = value;
}

void Foam::Vector_par::checkCompatible(const Vector_par& vec) const
{
    if(this->getLocalSize()!=vec.getLocalSize())
        FatalErrorInFunction<<"Vector local size mismatch"<<exit(FatalError);
    if(this->getGlobalSize()!=vec.getGlobalSize())
        FatalErrorInFunction<<"Vector global size mismatch"<<exit(FatalError);
}

Foam::Vector_par Foam::LinearSolver_par::solve
(
    const Vector_par& b
)
{
    Vector_par x(b);
    x.fill(0);
    return solve(b,x);
}

Foam::Vector_par Foam::LinearSolver_par::solve
(
    const Vector_par& b,
    const Vector_par& x
)
{
    this->b = b;
    this->x = x;
    initIteration();
    while(!step() || iteration<minIteration)
    {}
    return this->x;
}

void Foam::GaussSeidel::initIteration()
{}

bool Foam::GaussSeidel::step()
{   
    iteration++;
    
    Vector_par sum_x = A*x;
    for(label localRow=0; localRow<x.getLocalSize().second; localRow++)
    {
        scalar Akk_x = A.evalOnDiagonal(x,localRow);
        sum_x[localRow] -= Akk_x;
    }
    Vector_par AxGs = b-sum_x;
    for(label localRow=0; localRow<x.getLocalSize().second; localRow++)
    {
        scalar Akk = A.getDiagonal(localRow);
        AxGs[localRow] /= Akk;
    }
    x = AxGs;

    Vector_par r = A*x -b;
    if(r.norm2()<tolerance)
        return true;
    
    return false;
}

void Foam::ConjugateGradient::initIteration()
{
    r = b - A*x;
    d = r;
}

bool Foam::ConjugateGradient::step()
{   
    iteration++;
    
    Vector_par z = A*d;
    
    scalar alpha = (r&r)/(d&z);
    x = x + d*alpha;
    Vector_par r_prev = r;
    r = r - z*alpha;
    if(r.norm2()<tolerance)
        return true;
    
    scalar beta = (r&r)/(r_prev&r_prev);
    d = r + d*beta;
    return false;
}

bool Foam::BiCGSTAB::step()
{   
    iteration++;
    
    scalar rho = d&r;
    if(rho==0)
        FatalErrorInFunction<<"BiCGSTAB fail"<<exit(FatalError);
    
    Vector_par p;
    if(iteration==1)
        p = r;
    else
    {
        scalar beta = (rho/this->rho.back())*(this->alpha.back()/this->w.back());
        p = r + (this->p.back() - this->v.back() * this->w.back()) * beta;
    }
    
    this->v.push(A*p);
    this->alpha.push(rho/(d&this->v.back()));
    Vector_par s = r - this->v.back()*this->alpha.back();
    
    scalar norm_s = s.norm2();
    if(norm_s<tolerance)
    {
        x = x + p*this->alpha.back();
        return true;
    }
    
    Vector_par t = A*s;
    this->w.push((t&s)/(t&t));
    x = x + p*this->alpha.back() + s*this->w.back();
    r = s - t*this->w.back();
    
    if(r.norm2()<tolerance)
        return true;
    else
        if(this->w.back()==0)
            FatalErrorInFunction<<"Error: w=0"<<exit(FatalError);
    
    this->rho.push(rho);
    this->p.push(p);

    if(this->rho.size()>1)
        this->rho.pop();
    if(this->alpha.size()>1)
        this->alpha.pop();
    if(this->w.size()>1)
        this->w.pop();
    
    if(this->v.size()>1)
        this->v.pop();
    if(this->p.size()>1)
        this->p.pop();

    return false;
}

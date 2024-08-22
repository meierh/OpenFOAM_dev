#include "LinearSystem.H"

void Foam::compute_Ay
(
    const gismo::gsMatrix<scalar>& A,
    gismo::gsMatrix<scalar>& Ay,
    scalar y
)
{
    Ay = A;
    for(label row=0; row<Ay.rows(); row++)
    {
        for(label col=0; col<Ay.cols(); col++)
        {
            scalar& value = Ay(row,col);
            value = 0;
            if(row != col)
            {
                bool cond1 = std::abs(A(row,col))<y*std::min(std::abs(A(row,row)),std::abs(A(col,col)));
                bool cond2 = A(row,col)*A(row,row)>0;
                bool cond3 = (A(row,col)==A(col,row)) && (A(row,col)*A(col,col)>0);
                if(cond1 || cond2 || cond3)
                    value = A(row,col);
            }
        }
    }
}

void Foam::compute_ATy
(
    const gismo::gsMatrix<scalar>& A,
    gismo::gsMatrix<scalar>& ATy,
    scalar y
)
{
    gismo::gsMatrix<scalar> Ay;
    compute_Ay(A,Ay,y);
    
    gismo::gsMatrix<scalar> rs_Ay = Ay;
    for(label row=0; row<rs_Ay.rows(); row++)
    {
        scalar rowSum = 0;
        for(label col=0; col<Ay.cols(); col++)
        {
            rowSum += rs_Ay(row,col);
            rs_Ay(row,col) = 0;
        }
        rs_Ay(row,row) = rowSum;
    }
    
    ATy = A - Ay + rs_Ay;
}

void Foam::compute_c
(
    const gismo::gsMatrix<scalar>& A,
    scalar alpha,
    std::unordered_set<label>& c
)
{
    for(label row=0; row<A.rows(); row++)
        c.insert(row);
    
    for(label i=0; i<A.rows(); i++)
    {
        if(c.find(i)!=c.end())
        {
            for(label j=0; j<A.cols(); j++)
            {
                if(i!=j)
                {
                    if(std::abs(A(i,j))>=alpha*std::min(std::abs(A(i,i)),std::abs(A(j,j))))
                    {
                        c.erase(j);
                    }
                }
            }
        }
    }
}

void Foam::compute_fa
(
    const gismo::gsMatrix<scalar>& A,
    const std::unordered_set<label>& c,
    scalar alpha,
    scalar tau,
    std::unordered_set<label>& fa
)
{
    gismo::gsMatrix<scalar> ATalpha;
    compute_ATy(A,ATalpha,alpha);
    
    for(label i=0; i<A.rows(); i++)
    {
        if(c.find(i)==c.end())
        {
            scalar sum_Col_ATalpha = 0;
            for(label j : c)
            {
                sum_Col_ATalpha += ATalpha(i,j);
            }
            
            if(std::abs(sum_Col_ATalpha) >= tau*std::abs(A(i,i)))
            {
                fa.insert(i);
            }
        }
    }
}

void Foam::compute_fb
(
    const gismo::gsMatrix<scalar>& A,
    std::unordered_set<label>& c,
    const std::unordered_set<label>& fa,
    scalar alpha,
    scalar delta,
    std::unordered_set<label>& fb
)
{
    fb.clear();
    gismo::gsMatrix<scalar> ATalpha;
    compute_ATy(A,ATalpha,alpha); 

    for(label i=0; i<A.rows(); i++)
    {
        if(c.find(i)==c.end() && fa.find(i)==fa.end())
        {
            scalar sum_Col_ATalpha = 0;
            for(label j : c)
            {
                sum_Col_ATalpha += ATalpha(i,j);
            }
            for(label j : fa)
            {
                sum_Col_ATalpha += ATalpha(i,j);
            }
            if(std::abs(sum_Col_ATalpha) >= delta*std::abs(A(i,i)))
            {
                fb.insert(i);
            }
            else
            {
                c.insert(i);
            }
        }
    }
}

void Foam::compute_f
(
    const std::unordered_set<label>& fa,
    const std::unordered_set<label>& fb,
    std::unordered_set<label>& f
)
{
    for(label e : fa)
        f.insert(e);
    for(label e : fb)
        f.insert(e);
}

void Foam::coarsening
(
    const gismo::gsMatrix<scalar>& A,
    scalar alpha,
    scalar tau,
    scalar delta
)
{
    std::unordered_set<label> c;
    compute_c(A,alpha,c);
    
    std::unordered_set<label> fa;
    compute_fa(A,c,alpha,tau,fa);
    
    std::unordered_set<label> fb;
    compute_fb(A,c,fa,alpha,delta,fb);
    
    std::unordered_set<label> f;
    compute_f(fa,fb,f);
}

void Foam::algebraicMultigrid
(
    const gismo::gsMatrix<scalar>& A,
    gismo::gsMatrix<scalar>& x,
    const gismo::gsMatrix<scalar>& b,
    scalar tol
)
{
    gismo::gsMatrix<scalar> Ay = A;
    for(label row=0; row<Ay.rows(); row++)
    {
        for(label col=0; col<Ay.cols(); col++)
        {
            scalar& value = Ay(row,col);
            value = 0;
            
        }
    }
}
    
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
    
    CSR_Matrix_par A_csr(A.rows(),0,A.rows(),A.cols(),false);
    Vector_par x_csr(A.rows(),0,A.rows(),false);
    Vector_par b_csr(A.rows(),0,A.rows(),false);
    for(label row=0; row<A.rows(); row++)
    {
        List<scalar> oneRow(A.cols());
        for(label col=0; col<A.cols(); col++)
            oneRow[col] = A(row,col);
        A_csr.addRow(oneRow);
        x_csr[row] = x(row,0);
        b_csr[row] = b(row,0);
    }
    
    BiCGSTAB solver(A_csr);
    x_csr = solver.solve(b_csr,x_csr);
    Vector_par resid = A_csr*x_csr-b_csr;
    scalar norm2_resid = resid.norm2();

    if(norm2_resid>1e-8)
    {
        Pout<<"A:"<<A_csr.to_string()<<Foam::endl;
        Pout<<"x:"<<x_csr.to_string()<<Foam::endl;
        Pout<<"b:"<<b_csr.to_string()<<Foam::endl;
        Pout<<"A:"<<(A_csr*x_csr).to_string()<<Foam::endl;
        Pout<<"norm2_resid:"<<norm2_resid<<Foam::endl;

        FatalErrorInFunction<<"Failed computing weights"<<exit(FatalError);
    }
    
    for(label row=0; row<A.rows(); row++)
    {
        x_csr[row] = x(row,0) = x_csr[row];
    }
}

void Foam::luDecomposition
(
    const gismo::gsMatrix<scalar>& A,
    gismo::gsMatrix<scalar>& L,
    gismo::gsMatrix<scalar>& U
)
{
    if(A.cols()!=A.rows())
        FatalErrorInFunction<<"Matrix must be square"<<exit(FatalError);
    const label n = A.cols();
    
    L = gismo::gsMatrix<scalar>(n,n);
    U = gismo::gsMatrix<scalar>(n,n);
    for(label c=0; c<n; c++)
    {
        for(label r=0; r<n; r++)
        {
            L(r,c) = 0;
            U(r,c) = 0;
        }
    }
    
    for (label i=0; i<n; i++) 
    {
        // Upper Triangular
        for (label k=i; k<n; k++)
        {
            // Summation of L(i, j) * U(j, k)
            scalar sum = 0;
            for (label j=0; j<i; j++)
                sum += (L(i,j) * U(j,k));

            // Evaluating U(i, k)
            U(i,k) = A(i,k) - sum;
        }

        // Lower Triangular
        for (label k=i; k<n; k++) 
        {
            if (i == k)
                L(i,i) = 1; // Diagonal as 1
            else 
            {
                // Summation of L(k, j) * U(j, i)
                scalar sum = 0;
                for (label j=0; j<i; j++)
                    sum += (L(k,j) * U(j,i));

                // Evaluating L(k, i)
                if(U(i,i)==0)
                    FatalErrorInFunction<<"Divide by zero error incoming!"<<exit(FatalError);
                L(k,i) = (A(k,i) - sum) / U(i,i);
            }
        }
    }
}

void Foam::svd
(
    const gismo::gsMatrix<scalar>& A,
    gismo::gsMatrix<scalar>& U,
    gismo::gsMatrix<scalar>& S,
    gismo::gsMatrix<scalar>& V
)
{
    const std::size_t M = A.rows();
    const std::size_t N = A.cols();
    
    gsl_matrix* A_gsl = gsl_matrix_alloc(M,N);
    for(std::size_t row=0; row<M; row++)
    {
        for(std::size_t col=0; col<N; col++)
        {
            gsl_matrix_set(A_gsl,row,col,A(row,col));
        }
    }

    gsl_matrix* U_gsl = nullptr;
    gsl_matrix* V_gsl = gsl_matrix_alloc(N,N);
    gsl_vector* S_gsl = gsl_vector_alloc(N);
    gsl_vector* work = gsl_vector_alloc(N);

    int errCode = gsl_linalg_SV_decomp(A_gsl,V_gsl,S_gsl,work);
    if(errCode)
        FatalErrorInFunction<<"Error in SVD"<<exit(FatalError);
    U_gsl = A_gsl;

    auto gsl_to_gismo_matrix = [](gsl_matrix* gsl_mat, gismo::gsMatrix<scalar>& gismo_mat)
    {
        gismo_mat = gismo::gsMatrix<scalar>(gsl_mat->size1,gsl_mat->size2);
        for(std::size_t row=0; row<gsl_mat->size1; row++)
        {
            for(std::size_t col=0; col<gsl_mat->size2; col++)
            {
                gismo_mat(row,col) = gsl_matrix_get(gsl_mat,row,col);
            }
        }
    };
    
    auto gsl_to_gismo_vector = [](gsl_vector* gsl_vec, gismo::gsMatrix<scalar>& gismo_vec)
    {
        gismo_vec = gismo::gsMatrix<scalar>(gsl_vec->size,0);
        for(std::size_t row=0; row<gsl_vec->size; row++)
        {
            gismo_vec(row,1) = gsl_vector_get(gsl_vec,row);
        }
    };
    
    gsl_to_gismo_matrix(U_gsl,U);
    gsl_to_gismo_matrix(V_gsl,V);
    gsl_to_gismo_vector(S_gsl,S);

    gsl_matrix_free(U_gsl);
    gsl_matrix_free(V_gsl);
    gsl_vector_free(S_gsl);
    gsl_vector_free(work);
}

void Foam::eig
(
    const gismo::gsMatrix<scalar>& A,
    gismo::gsMatrix<scalar>& eigenvalues,
    gismo::gsMatrix<scalar>& eigenvectors
)
{
    if(A.rows()!=A.cols())
        FatalErrorInFunction<<"Eigenvalues only for square matrix"<<exit(FatalError);
    
    std::size_t N = A.rows();
    
    gsl_matrix* A_gsl = gsl_matrix_alloc(N,N);
    for(std::size_t row=0; row<N; row++)
    {
        for(std::size_t col=0; col<N; col++)
        {
            gsl_matrix_set(A_gsl,row,col,A(row,col));
        }
    }
    gsl_vector* eval = gsl_vector_alloc(N);
    gsl_matrix* evec = gsl_matrix_alloc(N,N);
    gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(N);
    
    gsl_eigen_symmv(A_gsl,eval,evec,w);

    auto gsl_to_gismo_matrix = [](gsl_matrix* gsl_mat, gismo::gsMatrix<scalar>& gismo_mat)
    {
        gismo_mat = gismo::gsMatrix<scalar>(gsl_mat->size1,gsl_mat->size2);
        for(std::size_t row=0; row<gsl_mat->size1; row++)
        {
            for(std::size_t col=0; col<gsl_mat->size2; col++)
            {
                gismo_mat(row,col) = gsl_matrix_get(gsl_mat,row,col);
            }
        }
    };
    
    auto gsl_to_gismo_vector = [](gsl_vector* gsl_vec, gismo::gsMatrix<scalar>& gismo_vec)
    {
        gismo_vec = gismo::gsMatrix<scalar>(gsl_vec->size,1);
        for(std::size_t row=0; row<gsl_vec->size; row++)
        {
            gismo_vec(row,0) = gsl_vector_get(gsl_vec,row);
        }
    };
    
    gsl_to_gismo_matrix(evec,eigenvectors);
    gsl_to_gismo_vector(eval,eigenvalues);
    
    
    gsl_matrix_free(A_gsl);
    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    gsl_eigen_symmv_free(w);
}

void Foam::eig
(
    const CSR_Matrix_par& mat,
    gismo::gsMatrix<scalar>& eigenvalues,
    gismo::gsMatrix<scalar>& eigenvectors
)
{
    gismo::gsMatrix<scalar> gsMat;
    CSR_Matrix_par_To_gismo_Local(mat,gsMat);
    eig(gsMat,eigenvalues,eigenvectors);
}

Foam::scalar Foam::determinant
(
    const gismo::gsMatrix<scalar>& A
)
{
    //std::cout<<"A:"<<std::endl<<A<<std::endl;
    gismo::gsMatrix<scalar> eval,evec;
    eig(A,eval,evec);
    //std::cout<<"eval:"<<std::endl<<eval<<std::endl;
    //std::cout<<"evec:"<<std::endl<<evec<<std::endl;
    
    scalar det=1;
    for(label row=0; row<eval.rows(); row++)
        det *= eval(row,0);
    return det;    
}

Foam::scalar Foam::determinant
(
    const CSR_Matrix_par& mat
)
{
    gismo::gsMatrix<scalar> gsMat;
    CSR_Matrix_par_To_gismo_Local(mat,gsMat);
    return determinant(gsMat);
}

Foam::scalar Foam::condition
(
    const gismo::gsMatrix<scalar>& A
)
{
    gismo::gsMatrix<scalar> eval,evec;
    eig(A,eval,evec);
    
    std::vector<scalar> evalList;
    for(label row=0; row<eval.rows(); row++)
        evalList.push_back(std::abs(eval(row,0)));
    scalar maxEval = *std::max_element(evalList.begin(),evalList.end());
    scalar minEval = *std::min_element(evalList.begin(),evalList.end());
    
    return maxEval/minEval;    
}

Foam::scalar Foam::condition
(
    const CSR_Matrix_par& mat
)
{
    gismo::gsMatrix<scalar> gsMat;
    CSR_Matrix_par_To_gismo_Local(mat,gsMat);
    return condition(gsMat);
}

void Foam::CSR_Matrix_par_To_gismo_Local
(
    const CSR_Matrix_par& mat,
    gismo::gsMatrix<scalar>& gsMat
)
{
    label globalRows = mat.getGlobalRows();
    label globalCols = mat.getGlobalCols();
    gsMat = gismo::gsMatrix<scalar>(globalRows,globalCols);
    
    for(label col=0; col<mat.getGlobalCols(); col++)
    {
        Vector_par colExtract(mat.getLocalRows(),mat.getLocalRowStart(),mat.getGlobalRows(),mat.getGlobal());
        colExtract.fill(0);
        std::pair<label,label> localRowInfo = colExtract.getLocalSize();
        label localRowStart = localRowInfo.first;
        label localRows = localRowInfo.second;
        if(col >= localRowStart)
        {
            label idx = col-localRowStart;
            if(idx < localRows)
            {
                colExtract[idx] = 1;
            }
        }
        
        const Vector_par columnVec = mat*colExtract;
        List<scalar> globalColumn;
        columnVec.collectGlobal(globalColumn);
        
        if(globalColumn.size()!=gsMat.rows())
            FatalErrorInFunction<<"Size mismatch"<<Foam::endl;
        
        for(label row=0; row<gsMat.rows(); row++)
            gsMat(row,col) = globalColumn[row];
    }
}

void Foam::invertableInfo
(
    const gismo::gsMatrix<scalar>& A,
    gismo::gsMatrix<scalar>& eigenvalues,
    gismo::gsMatrix<scalar>& eigenvectors,
    scalar& determinant,
    scalar& condition
)
{
    eig(A,eigenvalues,eigenvectors);
    
    determinant = 1;
    for(label row=0; row<eigenvalues.rows(); row++)
        determinant *= eigenvalues(row,0);
    
    std::vector<scalar> evalList;
    for(label row=0; row<eigenvalues.rows(); row++)
        evalList.push_back(std::abs(eigenvalues(row,0)));
    scalar maxEval = *std::max_element(evalList.begin(),evalList.end());
    scalar minEval = *std::min_element(evalList.begin(),evalList.end());
    condition = maxEval/minEval;    
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

Foam::Vector_par::Vector_par
(
    label localRows,
    label localRowStart,
    label globalRows,
    bool global
):
localRowStart(localRowStart),
localRows(localRows),
globalRows(globalRows),
global(global)
{
    V.setSize(localRows);
    globalCheck();
}
            
Foam::Vector_par::Vector_par
(
    const Vector_par& vec
):
localRowStart(vec.localRowStart),
localRows(vec.localRows),
globalRows(vec.globalRows),
global(vec.global)
{
    V = vec.V;
}

Foam::Vector_par::Vector_par
(
    bool global
):
localRowStart(-1),
localRows(-1),
globalRows(-1),
global(global)
{}

Foam::scalar& Foam::Vector_par::operator [](int idx)
{
    if(idx<0 || idx>=V.size())
    {
        Info<<"idx:"<<idx<<Foam::endl;
        FatalErrorInFunction<<"Out of range"<<exit(FatalError);
    }
    return V[idx];
}

Foam::scalar Foam::Vector_par::operator [](int idx) const
{
    if(idx<0 || idx>=V.size())
    {
        Info<<"idx:"<<idx<<Foam::endl;
        FatalErrorInFunction<<"Out of range"<<exit(FatalError);
    }
    return V[idx];
}

Foam::scalar Foam::Vector_par::operator&
(
    const Vector_par& vec
) const
{
    //std::cout<<"vec:"<<std::endl<<vec.to_string()<<std::endl;
    checkCompatible(vec);
    
    if(global)
    {
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
    else
    {
        scalar localValue = 0;
        for(label i=0; i<localRows; i++)
            localValue += (*this)[i]*vec[i];
        return localValue;
    }
}

Foam::Vector_par Foam::Vector_par::operator*
(
    scalar scale
) const
{
    Vector_par result(*this);
    for(label i=0; i<localRows; i++)
        result[i] = (*this)[i]*scale;
    return result;
}

Foam::Vector_par Foam::Vector_par::operator+
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

Foam::Vector_par Foam::Vector_par::operator-
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

Foam::Vector_par& Foam::Vector_par::operator=
(
    const Vector_par& vec
)
{
    this->V = vec.V;
    this->localRowStart = vec.localRowStart;
    this->localRows = vec.localRows;
    this->globalRows = vec.globalRows;
    this->global = vec.global;

    return *this;
}

Foam::scalar Foam::Vector_par::norm2() const
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

void Foam::Vector_par::fill
(
    List<scalar> values
)
{
    if(values.size()!=localRows)
        FatalErrorInFunction<<"Mismatch"<<exit(FatalError);
    
    for(label i=0; i<localRows; i++)
        (*this)[i] = values[i];
}

void Foam::Vector_par::collectGlobal
(
    List<scalar>& globalVector
) const
{
    if(!global)
        FatalErrorInFunction<<"Can not collect in non global vector"<<exit(FatalError);

    globalVector.setSize(globalRows);
    List<List<scalar>> globalData(Pstream::nProcs());
    globalData[Pstream::myProcNo()] = V;
    Pstream::gatherList(globalData);
    Pstream::scatterList(globalData);

    // Total size check
    label totalSize = 0;
    for(label proc=0; proc<Pstream::nProcs(); proc++)
        totalSize += globalData[proc].size();
    if(totalSize!=globalRows)
    {
        Pout<<"globalData:"<<globalData<<Foam::endl;
        Pout<<"Collected Size:"<<totalSize<<Foam::endl;
        Pout<<"Stored Size:"<<globalRows<<Foam::endl;
        FatalErrorInFunction<<"Mismatch in total vector size"<<exit(FatalError);
    }
    
    int i=0;
    for(label proc=0; proc<Pstream::myProcNo(); proc++)
    {
        for(label localRow=0; localRow<globalData[proc].size(); localRow++)
        {
            globalVector[i] = globalData[proc][localRow];
            i++;
        }
    }
    if(i!=localRowStart)
        FatalErrorInFunction<<"Mismatch in local offset"<<exit(FatalError);
    for(label proc=Pstream::myProcNo(); proc<Pstream::nProcs(); proc++)
    {
        for(label localRow=0; localRow<globalData[proc].size(); localRow++)
        {
            globalVector[i] = globalData[proc][localRow];
            i++;
        }
    }
}

std::string Foam::Vector_par::to_string() const
{
    List<scalar> globalVector;
    if(global)
        collectGlobal(globalVector);
    else
        globalVector = V;
    
    std::string vector = to_metaDataString();
    for(scalar vec : globalVector)
    {
        vector.append(" "+std::to_string(vec)+"\n");
    }
    return vector;
}

std::string Foam::Vector_par::to_metaDataString() const
{
    std::string metaData = " ("+std::to_string(globalRows)+")\n";
    metaData.append("["+std::to_string(localRowStart)+"-"+std::to_string(localRows)+"]\n");    
    std::string boolStr = global?"true":"false";
    metaData.append("g:"+boolStr+"\n");
    return metaData;
}

void Foam::Vector_par::checkCompatible(const Vector_par& vec) const
{
    if(this->getLocalSize()!=vec.getLocalSize())
        FatalErrorInFunction<<"Vector local size mismatch"<<exit(FatalError);
    if(this->getGlobalSize()!=vec.getGlobalSize())
        FatalErrorInFunction<<"Vector global size mismatch"<<exit(FatalError);
    if(!(this->global==vec.global))
    {
        std::cout<<"this->global:"<<this->global<<std::endl;
        std::cout<<"vec.global:"<<vec.global<<std::endl;

        Info<<"this->global!=vec.global:"<<(this->global!=vec.global)<<Foam::endl;
        Info<<"this:"<<this->to_string()<<Foam::endl;
        Info<<"vec:"<<vec.to_string()<<Foam::endl;
        FatalErrorInFunction<<"Vector global local mismatch"<<exit(FatalError);
    }
}

void Foam::Vector_par::globalCheck() const
{
    if(global)
    {
        List<label> globalLocalRowStart(Pstream::nProcs());
        globalLocalRowStart[Pstream::myProcNo()] = localRowStart;
        Pstream::gatherList(globalLocalRowStart);
        
        List<label> globalLocalRows(Pstream::nProcs());
        globalLocalRows[Pstream::myProcNo()] = localRows;
        Pstream::gatherList(globalLocalRows);
        
        List<label> globalGlobalRows(Pstream::nProcs());
        globalGlobalRows[Pstream::myProcNo()] = globalRows;
        Pstream::gatherList(globalGlobalRows);
        
        if(Pstream::master())
        {
            for(label globalRows : globalGlobalRows)
                if(globalRows!=this->globalRows)
                    FatalErrorInFunction<<"Global row size mismatch"<<exit(FatalError);
            for(label proc=0; proc<globalLocalRowStart.size()-1; proc++)
            {
                if(globalLocalRowStart[proc]+globalLocalRows[proc]!=globalLocalRowStart[proc+1])
                    FatalErrorInFunction<<"Local row mismatch"<<exit(FatalError);
            }
            if(globalLocalRowStart.last()+globalLocalRows.last()!=this->globalRows)
                FatalErrorInFunction<<"Local rows to global rows mismatch"<<exit(FatalError);
        }
    }
    else
    {
        if(localRowStart!=0)
            FatalErrorInFunction<<"Local row start mismatch as local"<<exit(FatalError);
        if(localRows!=globalRows)
            FatalErrorInFunction<<"Local global row mismatch as local"<<exit(FatalError);
    }
}

Foam::CSR_Matrix_par::CSR_Matrix_par
(
    label localRows,
    label localRowStart,
    label globalRows,
    label globalCols,
    bool global
):
localRowStart(localRowStart),
localRows(localRows),
globalRows(globalRows),
globalCols(globalCols),
global(global),
writtenRow(0),
lastWrittenRowIndex(0)
{
    V.setSize(0);
    Col_Index.setSize(0);
    Row_Index.setSize(0);
    
    globalCheck();
    if(localRows<1)
        setComplete();
    else
        complete = false;
}

Foam::CSR_Matrix_par::CSR_Matrix_par
(
    const CSR_Matrix_par& mat
):
V(mat.V),
Col_Index(mat.Col_Index),
Row_Index(mat.Row_Index),
localRowStart(mat.localRowStart),
localRows(mat.localRows),
globalRows(mat.globalRows),
globalCols(mat.globalCols),
complete(mat.complete),
global(mat.global),
writtenRow(mat.writtenRow),
lastWrittenRowIndex(mat.lastWrittenRowIndex)
{}

Foam::CSR_Matrix_par::CSR_Matrix_par
(
    bool global
):
localRowStart(-1),
localRows(-1),
globalRows(-1),
globalCols(-1),
global(global),
writtenRow(0),
lastWrittenRowIndex(0)
{
    complete = true;
}

void Foam::CSR_Matrix_par::addRow(List<scalar> row)
{
    if(row.size()!=globalCols)
    {
        Pout<<"row.size():"<<row.size()<<Foam::endl;
        Pout<<"globalCols:"<<globalCols<<Foam::endl;
        Pout<<"row:"<<row<<Foam::endl;
        FatalErrorInFunction<<"Invalid row size"<<exit(FatalError);
    }
    
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
    if(complete)
        FatalErrorInFunction<<"Can not add to already complete matrix"<<exit(FatalError);

    if(writtenRow<localRows)
        writtenRow++;
    else
    {
        Info<<"writtenRow:"<<writtenRow<<"  localRows:"<<localRows<<Foam::endl;
        FatalErrorInFunction<<"Out of range row"<<exit(FatalError);
    }
    
    Row_Index.append(lastWrittenRowIndex);
    lastWrittenRowIndex+=row.size();
    for(label i=0; i<row.size(); i++)
    {
        std::pair<scalar,label>& entry = row[i];
        scalar value = entry.first;
        label col = entry.second;
        if(col>=globalCols || col<0)
        {
            Info<<"col:"<<col<<Foam::endl;
            Info<<"globalCols:"<<globalCols<<Foam::endl;
            FatalErrorInFunction<<"Col out of range"<<exit(FatalError);
        }
        V.append(value);
        Col_Index.append(col);
    }
    if(writtenRow==localRows)
        setComplete();
}

void Foam::CSR_Matrix_par::addRows(List<List<scalar>> rows)
{   
    if(rows.size()!=localRows)
        FatalErrorInFunction<<"Invalid row size"<<exit(FatalError);
    
    for(List<scalar> row : rows)
        addRow(row);
    //Info<<"Added rows to CSR_Matrix_par: ["<<localRows<<","<<localRowStart<<","<<globalRows<<","<<globalCols<<"]"<<Foam::endl;
}

Foam::Vector_par Foam::CSR_Matrix_par::operator*
(
    const Vector_par& vec
) const
{
    if(!complete || globalCols==-1)
        FatalErrorInFunction<<"Incomplete matrix"<<exit(FatalError);
    checkCompatible(vec);
    
    if(global)
    {
        List<scalar> globalVector;
        vec.collectGlobal(globalVector);
        
        Vector_par result(vec);
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
                value += V[ind]*globalVector[col];
            }
            result[k] = value;
        }
        return result;
    }
    else
    {
        Vector_par result(vec);
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
                value += V[ind]*vec[col];
            }
            result[k] = value;
        }
        return result;
    }
}

Foam::CSR_Matrix_par Foam::CSR_Matrix_par::operator*
(
    const CSR_Matrix_par& mat
) const
{
    FatalErrorInFunction<<"Not yet implemented"<<exit(FatalError);
    return CSR_Matrix_par();
}

Foam::scalar Foam::CSR_Matrix_par::getDiagonal
(
    label localRow
) const
{
    if(!complete || globalCols==-1)
        FatalErrorInFunction<<"Incomplete matrix"<<exit(FatalError);
    
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

Foam::scalar Foam::CSR_Matrix_par::evalOnDiagonal
(
    const Vector_par& vec,
    label localRow
) const
{
    if(!complete || globalCols==-1)
        FatalErrorInFunction<<"Incomplete matrix"<<exit(FatalError);
    checkCompatible(vec);
    
    scalar diag_A = getDiagonal(localRow);
    if(diag_A==0)
        return 0;
    else
        return diag_A*vec[localRow];
}

bool Foam::CSR_Matrix_par::isSquare() const
{
    return (globalRows==globalCols);// && globalRows>0;
}

Foam::CSR_Matrix_par Foam::CSR_Matrix_par::diagonalMatrix() const
{
    if(!isSquare())
        FatalErrorInFunction<<"Must be square!"<<exit(FatalError);
    
    CSR_Matrix_par result(localRows,localRowStart,globalRows,globalCols,global);
    for(label localRow=0; localRow<localRows; localRow++)
    {
        label rowStartInd = Row_Index[localRow];
        label rowEndInd;
        if(localRow<localRows-1)
            rowEndInd = Row_Index[localRow+1];
        else
            rowEndInd = V.size();
        
        DynamicList<std::pair<scalar,label>> oneRow;        
        for(label Vind=rowStartInd; Vind<rowEndInd; Vind++)
        {
            label col = Col_Index[Vind];
            if(col==localRow+localRowStart)
                oneRow.append({V[Vind],col});
        }
        if(oneRow.size()>1)
            FatalErrorInFunction<<"Diagonal can not be larger than one"<<exit(FatalError);
        result.addRow(oneRow);
    }
    return result;
}

Foam::CSR_Matrix_par Foam::CSR_Matrix_par::offDiagonalMatrix() const
{
    CSR_Matrix_par result(localRows,localRowStart,globalRows,globalCols,global);
    for(label localRow=0; localRow<localRows; localRow++)
    {
        label rowStartInd = Row_Index[localRow];
        label rowEndInd;
        if(localRow<localRows-1)
            rowEndInd = Row_Index[localRow+1];
        else
            rowEndInd = V.size();
        
        DynamicList<std::pair<scalar,label>> oneRow;        
        for(label Vind=rowStartInd; Vind<rowEndInd; Vind++)
        {
            label col = Col_Index[Vind];
            if(col!=localRow+localRowStart)
                oneRow.append({V[Vind],col});
        }
        result.addRow(oneRow);
    }
    return result;
}

std::string Foam::CSR_Matrix_par::to_string() const
{
    label processCount;
    if(global)
        processCount = Pstream::nProcs();
    else
        processCount = 1;
    
    List<DynamicList<scalar>> global_V(processCount);
    List<DynamicList<label>> global_Col_Index(processCount);
    List<DynamicList<label>> global_Row_Index(processCount);
    
    if(global)
    {   
        global_V[Pstream::myProcNo()] = V;
        global_Col_Index[Pstream::myProcNo()] = Col_Index;
        global_Row_Index[Pstream::myProcNo()] = Row_Index;
    
        Pstream::gatherList(global_V);
        Pstream::gatherList(global_Col_Index);
        Pstream::gatherList(global_Row_Index);

        Pstream::scatterList(global_V);
        Pstream::scatterList(global_Col_Index);
        Pstream::scatterList(global_Row_Index);
    }
    else
    {
        global_V[0] = V;
        global_Col_Index[0] = Col_Index;
        global_Row_Index[0] = Row_Index;
    }
                
    std::string matrix = to_metaDataString();

    for(label proc=0; proc<processCount; proc++)
    {
        DynamicList<scalar>& procV = global_V[proc];
        DynamicList<label>& procCol_Index = global_Col_Index[proc];
        DynamicList<label>& procRow_Index = global_Row_Index[proc];
        /*
        if(Pstream::master())
        {
            Pout<<proc<<"  procV:"<<procV<<Foam::endl;
            Pout<<proc<<"  procCol_Index:"<<procCol_Index<<Foam::endl;
            Pout<<proc<<"  procRow_Index:"<<procRow_Index<<Foam::endl;
        }
        */
        
        for(label localRow=0; localRow<procRow_Index.size(); localRow++)
        {
            /*
            if(Pstream::master())
            {
                Pout<<proc<<"    localRow:"<<localRow<<Foam::endl;
                Pout<<proc<<"    localRows:"<<localRows<<Foam::endl;
            }
            */
            
            label rowStartInd = procRow_Index[localRow];
            label rowEndInd;
            if(localRow<procRow_Index.size()-1)
                rowEndInd = procRow_Index[localRow+1];
            else
                rowEndInd = procV.size();
            
            /*
            if(Pstream::master())
            {
                Pout<<proc<<"    rowStartInd:"<<rowStartInd<<Foam::endl;
                Pout<<proc<<"    --rowEndInd:"<<rowEndInd<<Foam::endl;
                Pout<<proc<<"    globalCols:"<<globalCols<<Foam::endl;
            }
            */
            
            List<scalar> fullRow(globalCols);
            for(label c=0; c<globalCols; c++)
                fullRow[c] = 0;
            
            for(label Vind=rowStartInd; Vind<rowEndInd; Vind++)
            {
                label col = procCol_Index[Vind];
                fullRow[col] = procV[Vind];
            }

            /*
            if(Pstream::master())
            {
                Pout<<proc<<"    fullRow:"<<fullRow<<Foam::endl;
            }
            */
        
            for(label c=0; c<globalCols; c++)
            {
                char numString[10];
                std::sprintf(numString,"%3.2e",fullRow[c]);
                matrix.append(" ");
                matrix.append(numString);
            }
            matrix.append("\n");

            /*
            if(Pstream::master())
            {
                Pout<<proc<<"    matrix:"<<matrix<<Foam::endl;
            }
            */
        }
        /*
        if(Pstream::master())
            Pout<<"DONE ------------"<<proc<<Foam::endl;
        */
    }    
    //Pout<<"------------ DONE ------------"<<Foam::endl;
    return matrix;
}

std::string Foam::CSR_Matrix_par::to_metaDataString() const
{
    std::string metaData = " ("+std::to_string(globalRows)+","+std::to_string(globalCols)+")\n";
    metaData.append("["+std::to_string(localRowStart)+"-"+std::to_string(localRows)+"]\n");    
    std::string boolStr = global?"true":"false";
    metaData.append("g:"+boolStr+"\n");
    return metaData;
}

void Foam::CSR_Matrix_par::checkCompatible
(
    const Vector_par& vec
) const
{
    if(globalCols!=vec.getGlobalSize())
    {
        Pout<<"globalCols:"<<globalCols<<Foam::endl;
        Pout<<"vec.getGlobalSize():"<<vec.getGlobalSize()<<Foam::endl;
        FatalErrorInFunction<<"Matrix Vector mismatch"<<exit(FatalError);
    }
    if(global!=vec.getGlobal())
    {
        Pout<<"this->global:"<<this->global<<Foam::endl;
        Pout<<"vec:"<<vec.getGlobal()<<Foam::endl;
        /*
        Pout<<"this:"<<this->to_string()<<Foam::endl;
        Pout<<"vec:"<<vec.to_string()<<Foam::endl;
        */
        FatalErrorInFunction<<"Matrix Vector global local mismatch"<<exit(FatalError);
    }
}

void Foam::CSR_Matrix_par::checkCompatible
(
    const CSR_Matrix_par& mat
) const
{
    if(globalCols!=mat.globalRows)
        FatalErrorInFunction<<"Matrix matrix mismatch"<<exit(FatalError);
    if(this->global!=mat.global)
        FatalErrorInFunction<<"Matrix matrix mismatch / global local"<<exit(FatalError);
}

void Foam::CSR_Matrix_par::globalCheck() const
{
    if(global)
    {
        List<label> globalLocalRowStart(Pstream::nProcs());
        globalLocalRowStart[Pstream::myProcNo()] = localRowStart;
        Pstream::gatherList(globalLocalRowStart);
        
        List<label> globalLocalRows(Pstream::nProcs());
        globalLocalRows[Pstream::myProcNo()] = localRows;
        Pstream::gatherList(globalLocalRows);
        
        List<label> globalGlobalRows(Pstream::nProcs());
        globalGlobalRows[Pstream::myProcNo()] = globalRows;
        Pstream::gatherList(globalGlobalRows);
        
        List<label> globalGlobalCols(Pstream::nProcs());
        globalGlobalCols[Pstream::myProcNo()] = globalCols;
        Pstream::gatherList(globalGlobalCols);
        
        if(Pstream::master())
        {
            for(label globalRows : globalGlobalRows)
                if(globalRows!=this->globalRows)
                    FatalErrorInFunction<<"Global row size mismatch"<<exit(FatalError);
            for(label globalCols : globalGlobalCols)
                if(globalCols!=this->globalCols)
                    FatalErrorInFunction<<"Global col size mismatch"<<exit(FatalError);
            for(label proc=0; proc<globalLocalRowStart.size()-1; proc++)
            {
                if(globalLocalRowStart[proc]+globalLocalRows[proc]!=globalLocalRowStart[proc+1])
                    FatalErrorInFunction<<"Local row mismatch"<<exit(FatalError);
            }
            if(globalLocalRowStart.last()+globalLocalRows.last()!=this->globalRows)
            {
                Pout<<"globalLocalRowStart.last():"<<globalLocalRowStart.last()<<Foam::endl;
                Pout<<"globalLocalRows.last():"<<globalLocalRows.last()<<Foam::endl;
                Pout<<"this->globalRows:"<<this->globalRows<<Foam::endl;
                FatalErrorInFunction<<"Local rows to global rows mismatch"<<exit(FatalError);
            }
        }
    }
    else
    {
        if(localRowStart!=0)
            FatalErrorInFunction<<"Local row start mismatch as local"<<exit(FatalError);
        if(localRows!=globalRows)
            FatalErrorInFunction<<"Local global row mismatch as local"<<exit(FatalError);
    }
}

void Foam::CSR_Matrix_par::setComplete()
{
    complete = true;
    if(Row_Index.size()!=localRows)
        FatalErrorInFunction<<"Local row to Row_Index mismatch:"<<localRows<<"!="<<Row_Index.size()<<exit(FatalError);
    if(V.size()!=Col_Index.size())
        FatalErrorInFunction<<"V and Col_Index mismatch"<<exit(FatalError);
    globalCheck();
}

Foam::CSR_DiagMatrix_par::CSR_DiagMatrix_par
(
    label localSize,
    label localStart,
    label globalSize,
    bool global
):
CSR_Matrix_par(localSize,localStart,globalSize,globalSize,global)
{}

Foam::CSR_DiagMatrix_par::CSR_DiagMatrix_par
(
    bool global
):
CSR_Matrix_par(global)
{}

Foam::CSR_DiagMatrix_par::CSR_DiagMatrix_par
(
    const Vector_par& vec
):
CSR_Matrix_par(vec.getLocalSize().second,vec.getLocalSize().first,vec.getGlobalSize(),vec.getGlobalSize(),vec.getGlobal())
{
    for(label localRow=0; localRow<localRows; localRow++)
    {
        addRow(vec[localRow]);
    }
}

void Foam::CSR_DiagMatrix_par::addRow(scalar element)
{
    CSR_Matrix_par::addRow(List<std::pair<scalar,label>>(1,{element,writtenRow}));
}

void Foam::CSR_DiagMatrix_par::addRows(List<scalar> diagonal)
{
    if(diagonal.size()!=localRows)
        FatalErrorInFunction<<"Invalid row size"<<exit(FatalError);
    for(scalar element : diagonal)
        addRow(element);
}

Foam::CSR_Matrix_par Foam::CSR_DiagMatrix_par::operator*
(
    const CSR_Matrix_par& mat
) const
{
    if(!complete || globalCols==-1)
        FatalErrorInFunction<<"Incomplete matrix"<<exit(FatalError);
    checkCompatible(mat);

    CSR_Matrix_par result(mat);
    
    label rowStart = 0;
    label rowEnd = 0;
    for(label r=0; r<result.Row_Index.size(); r++)
    {
        rowStart = result.Row_Index[r];
        rowEnd = ((r+1)<result.Row_Index.size()) ? result.Row_Index[r+1] : result.V.size();
        
        for(label c=rowStart; c<rowEnd; c++)
        {
            result.V[c] *= V[r];
        }
    }
    
    return result;
}

Foam::Vector_par Foam::CSR_DiagMatrix_par::operator*
(
    const Vector_par& vec
) const
{
    if(!complete || globalCols==-1)
        FatalErrorInFunction<<"Incomplete matrix"<<exit(FatalError);
    checkCompatible(vec);
            
    Vector_par result(vec);
    for(label k=0; k<this->localRows; k++)
        result[k] *= V[k];
    return result;
}

Foam::scalar Foam::CSR_DiagMatrix_par::getDiagonal
(
    label localRow
) const
{
    if(!complete || globalCols==-1)
        FatalErrorInFunction<<"Incomplete matrix"<<exit(FatalError);
    
    return V[localRow];
}

Foam::Vector_par Foam::LinearSolver_par::solve
(
    const Vector_par& b
)
{
    A.checkCompatible(b);
    
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
    A.checkCompatible(b);
    A.checkCompatible(x);
    b.checkCompatible(x);
    
    this->b = b;
    this->x = x;
    
    if(x.getGlobalSize()==0)
        return this->x;
    
    initIteration();
    while(!step() || iteration<minIteration)
    {
        if(iteration > 1e7)
            Info<<"High iteration number in linear solver"<<Foam::endl;
        if(norm2_r > 1e100)
            throw std::runtime_error("Linear solver diverged");
    }
    return this->x;
}

void Foam::Jacobi::initIteration()
{
    OffDiagonal = A.offDiagonalMatrix();
}

bool Foam::Jacobi::step()
{
    iteration++;
    
    Vector_par sum_x = OffDiagonal*x;
    Vector_par AxGs = b-sum_x;
    
    for(label localRow=0; localRow<x.getLocalSize().second; localRow++)
    {
        scalar Akk = A.getDiagonal(localRow);
        AxGs[localRow] /= Akk;
    }
    x = AxGs;

    r = A*x -b;
    norm2_r = r.norm2();
    if(norm2_r<tolerance)
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
    norm2_r = r.norm2();
    if(norm2_r<tolerance)
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
    
    norm2_r = r.norm2();
    if(norm2_r<tolerance)
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

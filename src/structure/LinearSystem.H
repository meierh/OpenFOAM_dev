#ifndef LINEARSYSTEM_H
#define LINEARSYSTEM_H

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include "activeRodMesh.h"
#pragma GCC diagnostic pop 

#include "gsl/gsl_linalg.h"
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_cblas.h"

#include "fvMesh.H"
#include <limits>
#include <queue>
#include <unordered_set>

namespace Foam
{
    void compute_Ay
    (
        const gismo::gsMatrix<scalar>& A,
        gismo::gsMatrix<scalar>& Ay,
        double y
    );
    
    void compute_ATy
    (
        const gismo::gsMatrix<scalar>& A,
        gismo::gsMatrix<scalar>& ATy,
        double y
    );
    
    void compute_c
    (
        const gismo::gsMatrix<scalar>& A,
        scalar alpha,
        std::unordered_set<label>& c
    );

    void compute_fa
    (
        const gismo::gsMatrix<scalar>& A,
        const std::unordered_set<label>& c,
        scalar alpha,
        scalar tau,
        std::unordered_set<label>& fa
    );

    void compute_fb
    (
        const gismo::gsMatrix<scalar>& A,
        std::unordered_set<label>& c,
        const std::unordered_set<label>& fa,
        scalar alpha,
        scalar delta,
        std::unordered_set<label>& fb
    );

    void compute_f
    (
        const std::unordered_set<label>& fa,
        const std::unordered_set<label>& fb,
        std::unordered_set<label>& f
    );

    void coarsening
    (
        const gismo::gsMatrix<scalar>& A,
        scalar alpha,
        scalar tau,
        scalar delta
    );
    
    void algebraicMultigrid
    (
        const gismo::gsMatrix<scalar>& A,
        gismo::gsMatrix<scalar>& x,
        const gismo::gsMatrix<scalar>& b,
        scalar tol
    );
    
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
    
    void luDecomposition
    (
        const gismo::gsMatrix<scalar>& A,
        gismo::gsMatrix<scalar>& L,
        gismo::gsMatrix<scalar>& U
    );
    
    void svd
    (
        const gismo::gsMatrix<scalar>& A,
        gismo::gsMatrix<scalar>& U,
        gismo::gsMatrix<scalar>& S,
        gismo::gsMatrix<scalar>& V
    );
    
    void eig
    (
        const gismo::gsMatrix<scalar>& A,
        gismo::gsMatrix<scalar>& eigenvalues,
        gismo::gsMatrix<scalar>& eigenvectors
    );
    
    class CSR_Matrix_par;
    
    void eig
    (
        const CSR_Matrix_par& mat,
        gismo::gsMatrix<scalar>& eigenvalues,
        gismo::gsMatrix<scalar>& eigenvectors
    );
    
    scalar determinant
    (
        const gismo::gsMatrix<scalar>& A
    );
    
    scalar determinant
    (
        const CSR_Matrix_par& mat
    );
    
    scalar condition 
    (
        const gismo::gsMatrix<scalar>& A
    );
    
    scalar condition
    (
        const CSR_Matrix_par& mat
    );
    
    void CSR_Matrix_par_To_gismo_Local
    (
        const CSR_Matrix_par& mat,
        gismo::gsMatrix<scalar>& gsMat
    );
    
    void invertableInfo
    (
        const gismo::gsMatrix<scalar>& A,
        gismo::gsMatrix<scalar>& eigenvalues,
        gismo::gsMatrix<scalar>& eigenvectors,
        scalar& determinant,
        scalar& condition
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
    
    void linearSolve_BiCGSTAB_par
    (
        const List<List<scalar>>& K_inv,
        const List<List<scalar>>& A,
        List<scalar>& x,
        const List<scalar>& b,
        List<scalar>& x_initial,
        scalar tol=1e-10
    );
    
    class Vector_par
    {
        public:
            Vector_par(label localRows, label localRowStart, label globalRows, bool global=true);
            Vector_par(const Vector_par& vec);
            Vector_par(bool global=true);
            
            scalar& operator [](int idx);
            scalar operator [](int idx) const;
            scalar operator& (const Vector_par&) const;
            Vector_par operator* (scalar) const;
            Vector_par operator+ (const Vector_par&) const;
            Vector_par operator- (const Vector_par&) const;
            Vector_par& operator= (const Vector_par&);
            scalar norm2() const;
            void scale(scalar scale);
            
            std::pair<label,label> getLocalSize() const {return {localRowStart,localRows};}
            label getGlobalSize() const {return globalRows;}
            bool getGlobal() const {return global;}
            
            void fill(scalar value);
            void fill(List<scalar> values);
            
            void collectGlobal(List<scalar>& globalVector) const;
            
            std::string to_string() const;
            std::string to_metaDataString() const;
            
            void checkCompatible(const Vector_par& vec) const;
            void globalCheck() const;
            
        private:            
            List<scalar> V;
            label localRowStart;
            label localRows;
            label globalRows;
            bool global;
            bool useBlas = true;
    };
    
    class CSR_DiagMatrix_par;
    
    class CSR_Matrix_par
    {
        public:
            CSR_Matrix_par(label localRows, label localRowStart, label globalRows, label globalCols, bool global=true);
            CSR_Matrix_par(const CSR_Matrix_par& vec);
            CSR_Matrix_par(bool global=true);
            
            virtual void addRow(List<scalar> row);
            virtual void addRow(List<std::pair<scalar,label>> row);
            virtual void addRows(List<List<scalar>> rows);
            
            Vector_par operator* (const Vector_par& vec) const;
            CSR_Matrix_par operator* (const CSR_Matrix_par& mat) const;
            CSR_Matrix_par& operator= (const CSR_Matrix_par& mat) = default;
            scalar getDiagonal(label localRow) const;
            virtual scalar evalOnDiagonal(const Vector_par& vec, label localRow) const;
            bool isSquare() const;
            
            label getLocalRowStart() const {return localRowStart;}
            label getLocalRows() const {return localRows;}
            label getGlobalRows() const {return globalRows;}
            label getGlobalCols() const {return globalCols;}
            bool getComplete() const {return complete;}
            bool getGlobal() const {return global;}
            
            virtual CSR_Matrix_par diagonalMatrix() const;
            virtual CSR_Matrix_par offDiagonalMatrix() const;
            
            std::string to_string() const;
            std::string to_metaDataString() const;
            
            void checkCompatible(const Vector_par& vec) const;
            void checkCompatible(const CSR_Matrix_par& mat) const;
            void globalCheck() const;
            
            virtual ~CSR_Matrix_par() = default;
            
        protected:
            void setComplete();
            
            DynamicList<scalar> V;
            DynamicList<label> Col_Index;
            DynamicList<label> Row_Index;
            label localRowStart;
            label localRows;
            label globalRows;
            label globalCols;
            bool complete;
            bool global;
            
            label writtenRow;
            label lastWrittenRowIndex;
            
            bool useBlas = true;

            friend class CSR_DiagMatrix_par;
    };
    
    class CSR_DiagMatrix_par : public CSR_Matrix_par
    {
        public:
            CSR_DiagMatrix_par(label localSize, label localStart, label globalSize, bool global=true);
            CSR_DiagMatrix_par(bool global=true);
            CSR_DiagMatrix_par(const Vector_par& vec);
            
            void addRow(List<scalar> row){deleted();};
            void addRow(List<std::pair<scalar,label>> row){deleted();};
            void addRows(List<List<scalar>> rows){deleted();};
            void addRow(scalar element);
            void addRows(List<scalar> diagonal);
            
            CSR_Matrix_par operator* (const CSR_Matrix_par& mat) const;
            Vector_par operator* (const Vector_par& vec) const;
            scalar getDiagonal(label localRow) const;
            scalar evalOnDiagonal(const Vector_par& vec, label localRow) const {deleted(); return 0;};
            bool isSquare() const {return true;};
            
            CSR_Matrix_par diagonalMatrix() const {deleted(); return CSR_Matrix_par();};
            CSR_Matrix_par offDiagonalMatrix() const {deleted(); return CSR_Matrix_par();};
            
        private:
            void deleted() const {FatalErrorInFunction<<"Deleted!"<<exit(FatalError);}
    };
    
    class LinearSolver_par
    {
        public:
            LinearSolver_par
            (
                const CSR_Matrix_par& A,
                scalar tolerance=1e-10,
                label minIteration=0
            ):A(A),tolerance(tolerance),minIteration(minIteration)
            {}
            
            LinearSolver_par
            (
                const CSR_Matrix_par& A,
                CSR_Matrix_par K_inv,
                scalar tolerance=1e-10,
                label minIteration=0
            ):A(A),K_inv(K_inv),tolerance(tolerance),minIteration(minIteration)
            {}
            
            virtual ~LinearSolver_par() = default;

            
            Vector_par solve(const Vector_par& b);
            Vector_par solve(const Vector_par& b, const Vector_par& x);
            
        protected:
            virtual void initIteration() = 0;
            virtual bool step() = 0;
            
            const CSR_Matrix_par& A;
            const CSR_Matrix_par K_inv;
            scalar tolerance;
            label minIteration;
            label iteration=0;
            
            Vector_par x,b,r;
            scalar norm2_r=0;
    };

    class Jacobi : public LinearSolver_par
    {
        public:
            Jacobi
            (
                const CSR_Matrix_par& A, 
                scalar tolerance=1e-10,
                label minIteration=0
            ):LinearSolver_par(A,tolerance,minIteration)
            {}
                        
        protected:
            void initIteration() override;
            bool step() override;
            
            CSR_Matrix_par OffDiagonal;
    };
    
    class ConjugateGradient : public LinearSolver_par
    {
        public:
            ConjugateGradient
            (
                const CSR_Matrix_par& A, 
                scalar tolerance=1e-10,
                label minIteration=0
            ):LinearSolver_par(A,tolerance,minIteration)
            {}
                        
        protected:
            void initIteration() override;
            bool step() override;
            Vector_par d;
    };
    
    class BiCGSTAB : public ConjugateGradient
    {
        public:
            BiCGSTAB
            (
                const CSR_Matrix_par& A, 
                scalar tolerance=1e-10,
                label minIteration=0
            ):ConjugateGradient(A,tolerance,minIteration)
            {}
                        
        protected:
            bool step() override;
            std::queue<scalar> rho,alpha,w;
            std::queue<Vector_par> v,p;
    };
}

#endif

#ifndef LINEARSYSTEM_H
#define LINEARSYSTEM_H

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include "activeRodMesh.h"
#pragma GCC diagnostic pop 

#include "fvCFD.H"
#include <limits>
#include <queue>

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
            
            inline scalar& operator [](int idx){return V[idx];}
            inline scalar operator [](int idx) const {return V[idx];}
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
    };
    
    class CSR_Matrix_par
    {
        public:
            CSR_Matrix_par(label localRows, label localRowStart, label globalRows, label globalCols, bool global=true);
            CSR_Matrix_par(bool global=true);
            
            void addRow(List<scalar> row);
            void addRow(List<std::pair<scalar,label>> row);
            
            void addRows(List<List<scalar>> rows);
            
            Vector_par operator* (const Vector_par& vec) const;
            CSR_Matrix_par operator* (const CSR_Matrix_par& mat) const;
            scalar getDiagonal(label localRow) const;
            scalar evalOnDiagonal(const Vector_par& vec, label localRow) const;
            bool isSquare() const;
            
            CSR_Matrix_par diagonalMatrix() const;
            CSR_Matrix_par offDiagonalMatrix() const;
            
            std::string to_string() const;
            std::string to_metaDataString() const;
            
            void checkCompatible(const Vector_par& vec) const;
            void checkCompatible(const CSR_Matrix_par& vec) const;
            void globalCheck() const;
            
        private:
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
            
            Vector_par x,b;
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
            Vector_par r,d;
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

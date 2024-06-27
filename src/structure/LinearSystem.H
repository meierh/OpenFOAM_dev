#ifndef LINEARSYSTEM_H
#define LINEARSYSTEM_H

#include "activeRodMesh.h"
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
            Vector_par(label localRows, label localRowStart, label globalRows):
            localRows(localRows),localRowStart(localRowStart),globalRows(globalRows)
            {V = List<scalar>(localRows);}
            
            Vector_par(const Vector_par& vec):
            localRows(vec.localRows),localRowStart(vec.localRowStart),globalRows(vec.globalRows)
            {V = vec.V;}
            
            Vector_par():localRows(-1),localRowStart(-1),globalRows(-1){}
            
            inline scalar& operator [](int idx);
            inline scalar operator [](int idx) const;
            scalar operator& (const Vector_par&) const;
            Vector_par operator* (scalar) const;
            Vector_par operator+ (const Vector_par&) const;
            Vector_par operator- (const Vector_par&) const;
            scalar norm2() const;
            void scale(scalar scale);
            
            std::pair<label,label> getLocalSize() const {return {localRowStart,localRows};}
            label getGlobalSize() const {return globalRows;}
            
            void fill(scalar value);
            
        private:
            void checkCompatible(const Vector_par& vec) const;
            
            List<scalar> V;
            label localRowStart;
            label localRows;
            label globalRows;
    };
    
    class CSR_Matrix_par
    {
        public:
            CSR_Matrix_par(label localRows, label localRowStart, label globalRows, label globalCols):
            localRows(localRows),localRowStart(localRowStart),globalRows(globalRows),globalCols(globalCols),
            writtenRow(0),lastWrittenRowIndex(0){}
            CSR_Matrix_par(): localRows(-1),localRowStart(-1),globalRows(-1),globalCols(-1),writtenRow(0),lastWrittenRowIndex(0){}
            
            void addRow(List<scalar> row);
            void addRow(List<std::pair<scalar,label>> row);
            
            Vector_par operator* (const Vector_par& vec) const;
            CSR_Matrix_par operator* (const CSR_Matrix_par& mat) const;
            scalar getDiagonal(label localRow) const;
            scalar evalOnDiagonal(const Vector_par& vec, label localRow) const ;
            
        private:
            void checkCompatible(const Vector_par& vec) const;
            void checkCompatible(const CSR_Matrix_par& vec) const;
            
            DynamicList<scalar> V;
            DynamicList<label> Col_Index;
            DynamicList<label> Row_Index;
            label localRowStart;
            label localRows;
            label globalRows;
            label globalCols;
            
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
            
            Vector_par solve(const Vector_par& b);
            Vector_par solve(const Vector_par& b, const Vector_par& x);
            
        protected:
            virtual void initIteration();
            virtual bool step();
            
            const CSR_Matrix_par& A;
            const CSR_Matrix_par K_inv;
            scalar tolerance;
            label minIteration;
            label iteration=0;
            
            Vector_par x,b;
    };

    class GaussSeidel : public LinearSolver_par
    {
        public:
            GaussSeidel
            (
                const CSR_Matrix_par& A, 
                scalar tolerance=1e-10,
                label minIteration=0
            ):LinearSolver_par(A,tolerance,minIteration)
            {}
                        
        protected:
            virtual void initIteration();
            virtual bool step();
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
            virtual void initIteration();
            virtual bool step();
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
            virtual bool step();
            std::queue<scalar> rho,alpha,w;
            std::queue<Vector_par> v,p;
    };
}

#endif

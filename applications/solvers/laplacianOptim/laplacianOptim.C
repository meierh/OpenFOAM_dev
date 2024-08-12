/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    icoFoam

Description
    System
        Equation:           alpha * T'' = q
        BoundaryCondition:  T(0) = T_0
                            T'(1) = dT_1
    Transient solver for incompressible, laminar flow of Newtonian fluids.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "LinearSystem.H"
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    
    // ----------------- Domain ------------------
    Info<<"Domain"<<Foam::endl;
    scalar global_xStart = 0;
    scalar global_xEnd = 1;
    
    // ----------------- Numeric grid -----------------
    Info<<"Numeric grid"<<Foam::endl;
    label globalCellCount = 99;
    globalCellCount = ((globalCellCount/Pstream::nProcs())+1)*Pstream::nProcs();
    scalar h,cellWidth;
    h = cellWidth = (global_xEnd-global_xStart)/globalCellCount;
    label cellPerProcess = globalCellCount/Pstream::nProcs();
    
    label local_cellIndStart = cellPerProcess*Pstream::myProcNo();
    
    List<scalar> local_cellCentres(cellPerProcess);
    List<Pair<label>> rightNeighbour(cellPerProcess);
    List<Pair<label>> leftNeighbour(cellPerProcess);
    for(label loc_CellInd=0; loc_CellInd<cellPerProcess; loc_CellInd++)
    {
        if(loc_CellInd==0)
        {
            leftNeighbour[loc_CellInd] = {Pstream::myProcNo()-1,cellPerProcess-1};
            rightNeighbour[loc_CellInd] = {Pstream::myProcNo(),loc_CellInd+1};
        }
        else if(loc_CellInd==cellPerProcess-1)
        {
            leftNeighbour[loc_CellInd] = {Pstream::myProcNo(),loc_CellInd-1};
            rightNeighbour[loc_CellInd] = {Pstream::myProcNo()+1,0};
        }
        else
        {
            leftNeighbour[loc_CellInd] = {Pstream::myProcNo(),loc_CellInd-1};
            rightNeighbour[loc_CellInd] = {Pstream::myProcNo(),loc_CellInd+1};
        }
    }
    
    // -------------- Coefficients -----------------------
    Info<<"Coefficients"<<Foam::endl;
    List<List<scalar>> alpha(Pstream::nProcs());
    alpha[Pstream::myProcNo()].resize(cellPerProcess);
    List<List<scalar>> q(Pstream::nProcs());
    q[Pstream::myProcNo()].resize(cellPerProcess);
    for(label loc_CellInd=0; loc_CellInd<cellPerProcess; loc_CellInd++)
    {
        alpha[Pstream::myProcNo()][loc_CellInd] = 1;
        q[Pstream::myProcNo()][loc_CellInd] = 0;
    }
    Pstream::gatherList(alpha);
    Pstream::scatterList(alpha);
    Pstream::gatherList(q);
    Pstream::scatterList(q);
    
    //---------------- BoundaryCondition -----------------
    Info<<"BoundaryCondition"<<Foam::endl;
    scalar T_0 = 10;
    scalar dT_1 = -10;
    
    //----------------- Field -------------------
    Info<<"Field"<<Foam::endl;
    List<List<scalar>> T(Pstream::nProcs());
    T[Pstream::myProcNo()].resize(cellPerProcess,0);
    
    // ------- Create laplacian system -----------
    Info<<"Create laplacian system"<<Foam::endl;
    CSR_Matrix_par M(cellPerProcess,local_cellIndStart,globalCellCount,globalCellCount);
    Vector_par b(cellPerProcess,local_cellIndStart,globalCellCount);

    Info<<"M:"<<M.to_metaDataString()<<Foam::endl;
    
    for(label p=0; p<cellPerProcess; p++)
    {
        List<std::pair<scalar,label>> cellRow;
        scalar b_p = 0;
        if(Pstream::myProcNo()==0 && p==0)
        {
            cellRow = List<std::pair<scalar,label>>(2);
            scalar alpha_E = 0.5*(alpha[rightNeighbour[p].first()][rightNeighbour[p].second()]+alpha[Pstream::myProcNo()][p]);
            scalar alpha_P = alpha[Pstream::myProcNo()][p];
            cellRow[0] = {-((alpha_E/h)+2*(alpha_P/h)),p};
            cellRow[1] = {alpha_E/h,p+1};
            b_p = q[Pstream::myProcNo()][p]*cellWidth-(2*alpha_P)*T_0/h;
        }
        else if(Pstream::myProcNo()==Pstream::nProcs()-1 && p==cellPerProcess-1)
        {
            cellRow = List<std::pair<scalar,label>>(2);
            scalar alpha_W = 0.5*(alpha[leftNeighbour[p].first()][leftNeighbour[p].second()]+alpha[Pstream::myProcNo()][p]);
            scalar alpha_P = alpha[Pstream::myProcNo()][p];
            cellRow[0] = {alpha_W/h,p-1};
            cellRow[1] = {-(alpha_W/h),p};
            b_p = q[Pstream::myProcNo()][p]*cellWidth-alpha_P*dT_1;
        }
        else
        {
            cellRow = List<std::pair<scalar,label>>(3);
            scalar alpha_W = 0.5*(alpha[leftNeighbour[p].first()][leftNeighbour[p].second()]+alpha[Pstream::myProcNo()][p]);
            scalar alpha_E = 0.5*(alpha[rightNeighbour[p].first()][rightNeighbour[p].second()]+alpha[Pstream::myProcNo()][p]);
            cellRow[0] = {alpha_E/h,p-1};
            cellRow[1] = {-((alpha_E/h)+(alpha_W/h)),p};
            cellRow[2] = {alpha_W/h,p+1};
            b_p = q[Pstream::myProcNo()][p]*cellWidth;
        }
        M.addRow(cellRow);
        b[p] = b_p;
    }
    
    // ------------- Solver ------------------
    //Info<<"M:"<<M.to_string()<<Foam::endl;
    //Info<<"b:"<<b.to_string()<<Foam::endl;
    Info<<"Solver"<<Foam::endl;
    BiCGSTAB solver(M);
    Vector_par T_sol = solver.solve(b);
    
    // -------------- Remap ------------------
    Info<<"Remap"<<Foam::endl;
    for(label p=0; p<cellPerProcess; p++)
    {
        T[Pstream::myProcNo()][p] = T_sol[p];
    }
    Pstream::gatherList(T);
    Pstream::scatterList(T);

    Info<<T<<Foam::endl;    

    // ------- Create adjoint laplacian system -----------
    Info<<"Create adjoint system"<<Foam::endl;
    CSR_Matrix_par M_adj = M;
    Vector_par b_adj(b);
    const scalar T_aim = 10;
    std::function<scalar(scalar)> djdT = [&](scalar T)
    {
        scalar deltaT = T-T_aim;
        return deltaT*deltaT;
    };

    Info<<"M_adj:"<<M_adj.to_metaDataString()<<Foam::endl;
    
    for(label p=0; p<cellPerProcess; p++)
    {
        scalar b_adj_p = 0;
        if(Pstream::myProcNo()==0 && p==0)
        {
            b_adj_p = 0;
        }
        else if(Pstream::myProcNo()==Pstream::nProcs()-1 && p==cellPerProcess-1)
        {
            b_adj_p = - djdT(T[Pstream::myProcNo()][p]);
        }
        else
        {
            b_adj_p = 0;
        }
        b_adj[p] = b_adj_p;
    }
    
    // ------------- Solver ------------------
    //Info<<"M_adj:"<<M_adj.to_string()<<Foam::endl;
    //Info<<"b_adj:"<<b_adj.to_string()<<Foam::endl;
    Info<<"Adjoint Solver"<<Foam::endl;
    BiCGSTAB adjSolver(M_adj);
    Vector_par Adj_sol = adjSolver.solve(b_adj);
    
    return 0;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    OFtutorial11_modifyingTheMesh

Description

\*---------------------------------------------------------------------------*/
#include <math.h>
#include <memory>


#include <stdio.h>
#include <iostream>
//#include <direct.h>
#include <random>
#include <chrono>

/*
#include "activeRodMesh.h"
#include "rodOpt.h"
#include "rodTools.h"
#include "rodLatticePBC.h"
#include "rodImport.h"
#include "rodCScircle.h"
#include "rodCSrectangle.h"
#include "timing.h"
*/

#include "fvCFD.H"
#include "LinearSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // set up the case
    #include "setRootCase.H"

    // create the run time object
    Info<< "Create time\n" << endl;
    Time runTime
    (
        Time::controlDictName,
        args.rootPath(),
        args.caseName()
    );

    List<List<scalar>> NonSymmA1 = {{1,-1,1},{1,1,-2},{3,-4,-7}};
    List<scalar> NonSymmb1 = {-2,9,0};
    List<scalar> NonSymmx1 = {3,4,-1};
    
    List<List<scalar>> NonSymmA2 = {{1,1,1,0},{-2,-2,2,3.33},{1,-1,-1,0},{0,2.5,3,-1}};
    List<scalar> NonSymmb2 = {2,2,8,-11};
    List<scalar> NonSymmx2 = {5,-2,-1,3};
    
    List<List<scalar>> NonSymmA3 = {{0.7071,0, 0,-1,-0.866, 0, 0, 0,      0},
                                    {0.7071,0, 1, 0,   0.5, 0, 0, 0,      0},
                                    {0     ,1, 0, 0,     0,-1, 0, 0,      0},
                                    {0     ,0,-1, 0,     0, 0, 0, 0,      0},
                                    {0     ,0, 0, 0,     0, 0, 1, 0, 0.7071},
                                    {0     ,0, 0, 1,     0, 0, 0, 0,-0.7071},
                                    {0     ,0, 0, 0,0.8660, 1, 0,-1,      0},
                                    {0     ,0, 0, 0,  -0.5, 0,-1, 0,      0},
                                    {0     ,0, 0, 0,     0, 0, 0, 1, 0.7071}};
    List<scalar> NonSymmb3 = {1000,-1000,0,500,500,0,0,-500,0};
    
    List<List<scalar>> NonSymmA4 = {{5.  ,2, 1,  0,   0,  0,  0, 0,   0},
                                    {0.5 ,3, 2,0.1,   0,  0,  0, 0,   0},
                                    {0   ,2, 7,  3,   1,  0,  0, 0,   0},
                                    {0   ,0, 0,0.5,   0,  0,  0, 0,   0},
                                    {1   ,0, 0,  0, 3.1,  0,  0, 0, 0.4},
                                    {0   ,0, 0,  0,   0,1.2,0.3, 0,   0},
                                    {0   ,2, 0,  0,   0,  1,  9, 0,   1},
                                    {0.1 ,0, 0,  0,   0,  0,  5, 7,   1},
                                    {0   ,0, 3,0.5,   0,  0,  0, 1,  12}};
    List<scalar> NonSymmb4 = {64.66,16.39,71.66,7,17.1,0,-3.09,22.2,82};
    
    List<List<scalar>> SymmA5 = {{4,1},{1,3}};
    List<scalar> Symmb5 = {1,2};
    List<scalar> Symmx5 = {0.090909090909090909,0.636363636363636363636363};   
    
    List<List<scalar>> SymmA6 = {{ 4 ,-1, 0,  0, -1,  0},
                                 {-1 , 4,-1,  0,  0, -1},
                                 { 0 ,-1, 4, -1,  0,  0},
                                 { 0 , 0,-1,  4,  0,  0},
                                 {-1 , 0, 0,  0,  4, -1},
                                 { 0 ,-1, 0,  0, -1,  4}};
    List<scalar> Symmb6 = {1,1,1,1,1,1};
    
    List<List<scalar>> SymmA7 = {{ 4 ,-1,-1,  0},
                                 {-1 , 4, 0, -1},
                                 {-1 , 0, 4, -1},
                                 { 0 ,-1,-1,  4}};
    List<scalar> Symmb7 = {8,4,4,1};
    
    List<List<scalar>> SymmA8 = {{ 4 , -1,  0, -1,  0,  0},
                                 {-1 ,  4, -1,  0, -1,  0},
                                 { 0 , -1,  4,  0,  0, -1},
                                 {-1 ,  0,  0,  4, -1,  0},
                                 { 0 , -1,  0, -1,  4, -1},
                                 { 0 ,  0, -1,  0, -1,  4}};
    List<scalar> Symmb8 = {12,7,3,2,5,15};
    
    auto solutionTest = [](CSR_Matrix_par& A, Vector_par& b)
    {
        Info<<"Solution Test"<<Foam::endl;
        std::string A_str = A.to_string();
        std::string b_str = b.to_string();
        Info<<"Strings 1"<<Foam::endl;
        
        Jacobi solverJ(A);
        Vector_par xSolJacobi = solverJ.solve(b);
        std::string xSolJacobi_str = xSolJacobi.to_string();
        std::string AxSolJacobi_str = (A*xSolJacobi).to_string();
        Info<<"Strings 2"<<Foam::endl;
        
        ConjugateGradient solverCG(A);
        Vector_par xSolCG = solverCG.solve(b);
        std::string xSolCG_str = xSolCG.to_string();
        std::string AxSolCG_str = (A*xSolCG).to_string();
        Info<<"Strings 3"<<Foam::endl;
        
        BiCGSTAB solverBCG(A);
        Vector_par xSolBCG = solverBCG.solve(b);
        std::string xSolBCG_str = xSolBCG.to_string();
        std::string AxSolBCG_str = (A*xSolBCG).to_string();
        Info<<"Strings 4"<<Foam::endl;
        
        if(Pstream::master())
        {
            Info<<"A:"<<A_str<<Foam::endl;
            Info<<"b:"<<b_str<<Foam::endl;
            
            Info<<"xSolJacobi_str:"<<xSolJacobi_str<<Foam::endl;
            Info<<"AxSolJacobi_str:"<<AxSolJacobi_str<<Foam::endl;
            
            Info<<"xSolCG_str:"<<xSolCG_str<<Foam::endl;
            Info<<"AxSolCG_str:"<<AxSolCG_str<<Foam::endl;
            
            Info<<"xSolBCG_str:"<<xSolBCG_str<<Foam::endl;
            Info<<"AxSolBCG_str:"<<AxSolBCG_str<<Foam::endl;
        }
    };
    
    if(Pstream::parRun() && Pstream::nProcs()==2)
    {        
        CSR_Matrix_par A1,A2,A3,A4,A5,A6,A7,A8;
        Vector_par x1,x2,x3,x4,x5;
        Vector_par b1,b2,b3,b4,b5,b6,b7,b8;
        
        if(Pstream::myProcNo()==0)
        {
            A1 = CSR_Matrix_par(2,0,NonSymmA1.size(),NonSymmA1[0].size());
            A1.addRow(NonSymmA1[0]);
            A1.addRow(NonSymmA1[1]);
            x1 = Vector_par(2,0,NonSymmx1.size());
            x1[0] = NonSymmx1[0];
            x1[1] = NonSymmx1[1];
            b1 = Vector_par(2,0,NonSymmb1.size());
            b1[0] = NonSymmb1[0];
            b1[1] = NonSymmb1[1];
            
            A2 = CSR_Matrix_par(2,0,NonSymmA2.size(),NonSymmA2[0].size());
            A2.addRow(NonSymmA2[0]);
            A2.addRow(NonSymmA2[1]);
            x2 = Vector_par(2,0,NonSymmx2.size());
            x2[0] = NonSymmx2[0];
            x2[1] = NonSymmx2[1];
            b2 = Vector_par(2,0,NonSymmb2.size());
            b2[0] = NonSymmb2[0];
            b2[1] = NonSymmb2[1];
            
            A3 = CSR_Matrix_par(5,0,NonSymmA3.size(),NonSymmA3[0].size());
            A3.addRow(NonSymmA3[0]);
            A3.addRow(NonSymmA3[1]);
            A3.addRow(NonSymmA3[2]);
            A3.addRow(NonSymmA3[3]);
            A3.addRow(NonSymmA3[4]);
            b3 = Vector_par(5,0,NonSymmb3.size());
            b3[0] = NonSymmb3[0];
            b3[1] = NonSymmb3[1];
            b3[2] = NonSymmb3[2];
            b3[3] = NonSymmb3[3];
            b3[4] = NonSymmb3[4];
            
            A4 = CSR_Matrix_par(5,0,NonSymmA4.size(),NonSymmA4[0].size());
            A4.addRow(NonSymmA4[0]);
            A4.addRow(NonSymmA4[1]);
            A4.addRow(NonSymmA4[2]);
            A4.addRow(NonSymmA4[3]);
            A4.addRow(NonSymmA4[4]);
            b4 = Vector_par(5,0,NonSymmb4.size());
            b4[0] = NonSymmb4[0];
            b4[1] = NonSymmb4[1];
            b4[2] = NonSymmb4[2];
            b4[3] = NonSymmb4[3];
            b4[4] = NonSymmb4[4];
            
            A5 = CSR_Matrix_par(1,0,SymmA5.size(),SymmA5[0].size());
            A5.addRow(SymmA5[0]);
            x5 = Vector_par(1,0,Symmx5.size());
            x5[0] = Symmx5[0];
            b5 = Vector_par(1,0,Symmb5.size());
            b5[0] = Symmb5[0];
        
            A6 = CSR_Matrix_par(3,0,SymmA6.size(),SymmA6[0].size());
            A6.addRow(SymmA6[0]);
            A6.addRow(SymmA6[1]);
            A6.addRow(SymmA6[2]);
            b6 = Vector_par(3,0,Symmb6.size());
            b6[0] = Symmb6[0];
            b6[1] = Symmb6[1];
            b6[2] = Symmb6[2];  
        
            A7 = CSR_Matrix_par(2,0,SymmA7.size(),SymmA7[0].size());
            A7.addRow(SymmA7[0]);
            A7.addRow(SymmA7[1]);
            b7 = Vector_par(2,0,Symmb7.size());
            b7[0] = Symmb7[0];
            b7[1] = Symmb7[1];
        
            A8 = CSR_Matrix_par(3,0,SymmA8.size(),SymmA8[0].size());
            A8.addRow(SymmA8[0]);
            A8.addRow(SymmA8[1]);
            A8.addRow(SymmA8[2]);
            b8 = Vector_par(3,0,Symmb8.size());
            b8[0] = Symmb8[0];
            b8[1] = Symmb8[1];
            b8[2] = Symmb8[2];
        }
        else
        {
            A1 = CSR_Matrix_par(1,2,NonSymmA1.size(),NonSymmA1[0].size());
            A1.addRow(NonSymmA1[2]);
            x1 = Vector_par(1,2,NonSymmx1.size());
            x1[0] = NonSymmx1[2];
            b1 = Vector_par(1,2,NonSymmb1.size());
            b1[0] = NonSymmb1[2];
        
            A2 = CSR_Matrix_par(2,2,NonSymmA2.size(),NonSymmA2[0].size());
            A2.addRow(NonSymmA2[2]);
            A2.addRow(NonSymmA2[3]);
            x2 = Vector_par(2,2,NonSymmx2.size());
            x2[0] = NonSymmx2[2];
            x2[1] = NonSymmx2[3];
            b2 = Vector_par(2,2,NonSymmb2.size());
            b2[0] = NonSymmb2[2];
            b2[1] = NonSymmb2[3];
    
            A3 = CSR_Matrix_par(4,5,NonSymmA3.size(),NonSymmA3[0].size());
            A3.addRow(NonSymmA3[5]);
            A3.addRow(NonSymmA3[6]);
            A3.addRow(NonSymmA3[7]);
            A3.addRow(NonSymmA3[8]);
            b3 = Vector_par(4,5,NonSymmb3.size());
            b3[0] = NonSymmb3[5];
            b3[1] = NonSymmb3[6];
            b3[2] = NonSymmb3[7];
            b3[3] = NonSymmb3[8];
            
            A4 = CSR_Matrix_par(4,5,NonSymmA4.size(),NonSymmA4[0].size());
            A4.addRow(NonSymmA4[5]);
            A4.addRow(NonSymmA4[6]);
            A4.addRow(NonSymmA4[7]);
            A4.addRow(NonSymmA4[8]);
            b4 = Vector_par(4,5,NonSymmb4.size());
            b4[0] = NonSymmb4[5];
            b4[1] = NonSymmb4[6];
            b4[2] = NonSymmb4[7];
            b4[3] = NonSymmb4[8];
            
            A5 = CSR_Matrix_par(1,1,SymmA5.size(),SymmA5[0].size());
            A5.addRow(SymmA5[1]);
            x5 = Vector_par(1,1,Symmx5.size());
            x5[0] = Symmx5[1];
            b5 = Vector_par(1,1,Symmb5.size());
            b5[0] = Symmb5[1];
        
            A6 = CSR_Matrix_par(3,3,SymmA6.size(),SymmA6[0].size());
            A6.addRow(SymmA6[3]);
            A6.addRow(SymmA6[4]);
            A6.addRow(SymmA6[5]);
            b6 = Vector_par(3,3,Symmb6.size());
            b6[0] = Symmb6[3];
            b6[1] = Symmb6[4];
            b6[2] = Symmb6[5];  
        
            A7 = CSR_Matrix_par(2,2,SymmA7.size(),SymmA7[0].size());
            A7.addRow(SymmA7[2]);
            A7.addRow(SymmA7[3]);
            b7 = Vector_par(2,2,Symmb7.size());
            b7[0] = Symmb7[2];
            b7[1] = Symmb7[3];
        
            A8 = CSR_Matrix_par(3,3,SymmA8.size(),SymmA8[0].size());
            A8.addRow(SymmA8[3]);
            A8.addRow(SymmA8[4]);
            A8.addRow(SymmA8[5]);
            b8 = Vector_par(3,3,Symmb8.size());
            b8[0] = Symmb8[3];
            b8[1] = Symmb8[4];
            b8[2] = Symmb8[5];
        }
        
        // Matrix Vector multiplication
        Info<<"-------------------- Matrix vector multiplication ---------------------"<<Foam::endl;
        
        Vector_par A1b1 = A1*b1;
        std::string A1b1_str = A1b1.to_string();
        Vector_par A1x1 = A1*x1;
        std::string A1x1_str = A1x1.to_string();
        
        Vector_par A2b2 = A2*b2;
        std::string A2b2_str = A2b2.to_string();
        Vector_par A2x2 = A2*x2;
        std::string A2x2_str = A2x2.to_string();
        
        Vector_par A3b3 = A3*b3;
        std::string A3b3_str = A3b3.to_string();
        
        Vector_par A4b4 = A4*b4;
        std::string A4b4_str = A4b4.to_string();
        
        if(Pstream::master())
        {
            Info<<"A1b1_str:"<<A1b1_str<<Foam::endl;
            Info<<"A1x1_str:"<<A1x1_str<<Foam::endl;
            Info<<"A2b2_str:"<<A2b2_str<<Foam::endl;
            Info<<"A2x2_str:"<<A2x2_str<<Foam::endl;
            Info<<"A3b3_str:"<<A3b3_str<<Foam::endl;
            Info<<"A4b4_str:"<<A4b4_str<<Foam::endl;
        }

        //LGS Solution
        Info<<"----------------------------- LGS Solution -----------------------------"<<Foam::endl;
        auto solutionTest = [](CSR_Matrix_par& A, Vector_par& b)
        {
            std::string A_str = A.to_string();
            std::string b_str = b.to_string();
            
            Jacobi solverJ(A);
            Vector_par xSolJacobi = solverJ.solve(b);
            std::string xSolJacobi_str = xSolJacobi.to_string();
            std::string AxSolJacobi_str = (A*xSolJacobi).to_string();
            
            ConjugateGradient solverCG(A);
            Vector_par xSolCG = solverCG.solve(b);
            std::string xSolCG_str = xSolCG.to_string();
            std::string AxSolCG_str = (A*xSolCG).to_string();
            
            BiCGSTAB solverBCG(A);
            Vector_par xSolBCG = solverBCG.solve(b);
            std::string xSolBCG_str = xSolBCG.to_string();
            std::string AxSolBCG_str = (A*xSolBCG).to_string();
            
            if(Pstream::master())
            {
                Info<<"A:"<<A_str<<Foam::endl;
                Info<<"b:"<<b_str<<Foam::endl;
                
                Info<<"xSolJacobi_str:"<<xSolJacobi_str<<Foam::endl;
                Info<<"AxSolJacobi_str:"<<AxSolJacobi_str<<Foam::endl;
                
                Info<<"xSolCG_str:"<<xSolCG_str<<Foam::endl;
                Info<<"AxSolCG_str:"<<AxSolCG_str<<Foam::endl;
                
                Info<<"xSolBCG_str:"<<xSolBCG_str<<Foam::endl;
                Info<<"AxSolBCG_str:"<<AxSolBCG_str<<Foam::endl;
            }
        };
        
        solutionTest(A5,b5);
        solutionTest(A6,b6);
        solutionTest(A7,b7);
        solutionTest(A8,b8);
    }
    else if(Pstream::parRun() && Pstream::nProcs()==4)
    {        
        CSR_Matrix_par A8;
        Vector_par b8;

        if(Pstream::myProcNo()==0)
        {        
            A8 = CSR_Matrix_par(2,0,SymmA8.size(),SymmA8[0].size());
            A8.addRow(SymmA8[0]);
            A8.addRow(SymmA8[1]);
            b8 = Vector_par(2,0,Symmb8.size());
            b8[0] = Symmb8[0];
            b8[1] = Symmb8[1];
            
            CSR_Matrix_par A7(SymmA7.size(),0,SymmA7.size(),SymmA7[0].size(),false);
            A7.addRow(SymmA7[0]);
            A7.addRow(SymmA7[1]);
            A7.addRow(SymmA7[2]);
            A7.addRow(SymmA7[3]);
            Vector_par b7(Symmb7.size(),0,Symmb7.size(),false);
            b7[0] = Symmb7[0];
            b7[1] = Symmb7[1];
            b7[2] = Symmb7[2];
            b7[3] = Symmb7[3];
            
            Info<<"A7 local:"<<A7.to_string()<<Foam::endl;
            Info<<"b7 local:"<<b7.to_string()<<Foam::endl;
            
            Vector_par A7b7 = A7*b7;
            Info<<"A7b7 local:"<<A7b7.to_string()<<Foam::endl;
            solutionTest(A7,b7);
        }
        else if(Pstream::myProcNo()==1)
        {        
            A8 = CSR_Matrix_par(0,2,SymmA8.size(),SymmA8[0].size());
            b8 = Vector_par(0,2,Symmb8.size());
        }
        else if(Pstream::myProcNo()==2)
        {        
            A8 = CSR_Matrix_par(3,2,SymmA8.size(),SymmA8[0].size());
            A8.addRow(SymmA8[2]);
            A8.addRow(SymmA8[3]);
            A8.addRow(SymmA8[4]);
            b8 = Vector_par(3,2,Symmb8.size());
            b8[0] = Symmb8[2];
            b8[1] = Symmb8[3];
            b8[2] = Symmb8[4];
        }
        else
        {        
            A8 = CSR_Matrix_par(1,5,SymmA8.size(),SymmA8[0].size());
            A8.addRow(SymmA8[5]);
            b8 = Vector_par(1,5,Symmb8.size());
            b8[0] = Symmb8[5];
        }

        // Matrix Vector multiplication
        Info<<"-------------------- Matrix vector multiplication ---------------------"<<Foam::endl;
        Vector_par A8b8 = A8*b8;
        std::string A8b8_str = A8b8.to_string();
        if(Pstream::master())
        {
            Info<<"A8b8_str:"<<A8b8_str<<Foam::endl;
        }

        //LGS Solution
        Info<<"----------------------------- LGS Solution -----------------------------"<<Foam::endl;
        solutionTest(A8,b8);
    }
    else
    {
        CSR_Matrix_par A1(NonSymmA1.size(),0,NonSymmA1.size(),NonSymmA1[0].size());
        A1.addRows(NonSymmA1);
        Vector_par x1(NonSymmx1.size(),0,NonSymmx1.size());
        x1.fill(NonSymmx1);
        Vector_par b1(NonSymmb1.size(),0,NonSymmb1.size());
        b1.fill(NonSymmb1);
        
        CSR_Matrix_par A2(NonSymmA2.size(),0,NonSymmA2.size(),NonSymmA2[0].size());
        A2.addRows(NonSymmA2);
        Vector_par x2(NonSymmx2.size(),0,NonSymmx2.size());
        x2.fill(NonSymmx2);
        Vector_par b2(NonSymmb2.size(),0,NonSymmb2.size());
        b2.fill(NonSymmb2);
        
        CSR_Matrix_par A3(NonSymmA3.size(),0,NonSymmA3.size(),NonSymmA3[0].size());
        A3.addRows(NonSymmA3);
        Vector_par b3(NonSymmb3.size(),0,NonSymmb3.size());
        b3.fill(NonSymmb3);
        
        CSR_Matrix_par A4(NonSymmA4.size(),0,NonSymmA4.size(),NonSymmA4[0].size());
        A4.addRows(NonSymmA4);
        Vector_par b4(NonSymmb4.size(),0,NonSymmb4.size());
        b4.fill(NonSymmb4);
        
        CSR_Matrix_par A5(SymmA5.size(),0,SymmA5.size(),SymmA5[0].size());
        A5.addRows(SymmA5);
        Vector_par x5(Symmx5.size(),0,Symmx5.size());
        x5.fill(Symmx5);
        Vector_par b5(Symmb5.size(),0,Symmb5.size());
        b5.fill(Symmb5);
    
        CSR_Matrix_par A6(SymmA6.size(),0,SymmA6.size(),SymmA6[0].size());
        A6.addRows(SymmA6);
        Vector_par b6(Symmb6.size(),0,Symmb6.size());
        b6.fill(Symmb6);
        
        CSR_Matrix_par A7(SymmA7.size(),0,SymmA7.size(),SymmA7[0].size());
        A7.addRows(SymmA7);
        Vector_par b7(Symmb7.size(),0,Symmb7.size());
        b7.fill(Symmb7);
        
        CSR_Matrix_par A8(SymmA8.size(),0,SymmA8.size(),SymmA8[0].size());
        A8.addRows(SymmA8);
        Vector_par b8(Symmb8.size(),0,Symmb8.size());
        b8.fill(Symmb8);
        
        // Matrix Vector multiplication
        Info<<"-------------------- Matrix vector multiplication ---------------------"<<Foam::endl;
        Vector_par A1b1 = A1*b1;
        Info<<"A1b1:"<<A1b1.to_string()<<Foam::endl;
        Vector_par A1x1 = A1*x1;
        Info<<"A1x1:"<<A1x1.to_string()<<Foam::endl;
                
        Vector_par A2b2 = A2*b2;
        Info<<"A2b2:"<<A2b2.to_string()<<Foam::endl;
        Vector_par A2x2 = A2*x2;
        Info<<"A2x2:"<<A2x2.to_string()<<Foam::endl;

        Vector_par A3b3 = A3*b3;
        Info<<"A3b3:"<<A3b3.to_string()<<Foam::endl;
                
        Vector_par A4b4 = A4*b4;
        Info<<"A4b4:"<<A4b4.to_string()<<Foam::endl;
        
        Vector_par A8b8 = A8*b8;
        Info<<"A8b8_str:"<<A8b8.to_string();
                
        //LGS Solution
        Info<<"----------------------------- LGS Solution -----------------------------"<<Foam::endl;
        
        solutionTest(A5,b5);
        solutionTest(A6,b6);
        solutionTest(A7,b7);
        solutionTest(A8,b8);
    }

    Info << nl << "To best visualise the results, load the mesh and extract all patches" << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

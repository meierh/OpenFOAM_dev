/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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
    structureOptimTest

Description
    Testing of derivative rod methods

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "icoAdjointImmersedBoundary.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar RandomFloat(scalar lower, scalar upper)
{
    scalar random = static_cast<scalar>(rand()) / static_cast<scalar>(RAND_MAX);
    scalar diff = upper - lower;
    scalar delta = random * diff;
    return lower + delta;
}

label RandomInt(label lower, label upper)
{
    return static_cast<label>(RandomFloat(lower,upper));
}

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info<<"Start testing"<<Foam::endl;
    
    for(label run=0; run<1; run++)
    {
        scalar a0 = RandomFloat(0,1);
        label k = RandomInt(0,10);
        Info<<"----- Run:"<<run<<" -----:"<<k<<Foam::endl;

        std::vector<scalar> ak(k);
        std::vector<scalar> bk(k);
        for(label i=0; i<k; i++)
        {
            ak[i] = RandomFloat(0,1);
            bk[i] = RandomFloat(0,1);
        }
        scalar phase = RandomFloat(0,6.2);
        std::cout<<"a0:"<<a0<<std::endl;
        std::cout<<"ak: ";
        for(scalar a: ak)
            std::cout<<a<<" ";
        std::cout<<std::endl;
        std::cout<<"bk: ";
        for(scalar b: bk)
            std::cout<<b<<" ";
        std::cout<<std::endl;
        std::cout<<"phase:"<<phase<<std::endl;
        
        List<List<vector>> curveCoeffs(1);
        curveCoeffs[0] = List<vector>(4);
        scalar eps;
        eps = RandomFloat(-1,1);
        curveCoeffs[0][0] = vector(0+eps,0+eps,0+eps);
        eps = RandomFloat(-1,1);
        curveCoeffs[0][1] = vector(10+eps,10+eps,-1+eps);
        eps = RandomFloat(-1,1);
        curveCoeffs[0][2] = vector(20+eps,0+eps,1+eps);
        eps = RandomFloat(-1,1);
        curveCoeffs[0][3] = vector(40+eps,-10+eps,0+eps);
        Info<<"curveCoeffs:"<<curveCoeffs<<Foam::endl;
        
        std::vector<CrossSection> crossSecList = {CrossSection(a0,ak,bk,phase)};
        CrossSectionStructure testStructure(mesh,crossSecList,true);
        testStructure.setCurveCoeffs(curveCoeffs);
        testStructure.selfCheck();
    }
    
    Info<<"Done testing"<<Foam::endl;    
    return 0;
}



// ************************************************************************* //

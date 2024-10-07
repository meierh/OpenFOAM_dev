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

    solvers::icoAdjointImmersedBoundary solver(mesh,runTime);
    
    Info<<"Done testing"<<Foam::endl;    
    return 0;
}



// ************************************************************************* //

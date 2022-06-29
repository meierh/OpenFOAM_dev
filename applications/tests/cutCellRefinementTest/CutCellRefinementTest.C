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

//#include "dynamicRefineFvMesh.H"
#include "cutCellFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void UnitTests(int argc, char *argv[],Time& runTime)
{
    Info<<"---------------------------Unit Cut Refinement Test---------------------------"<<endl;
    UnitTest_cutCellFvMesh(argc,argv,runTime);
    Info<<"-------------------------Unit Test End-------------------------"<<endl;

}

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

    // disable post-processing etc.
    runTime.functionObjects().off();

    
    Foam::Info
    << "Create mesh for time = "
    << runTime.timeName() << Foam::nl << Foam::endl;

    UnitTests(argc,argv,runTime);
    
    Foam::Info<<"Test Nurbs Curve"<<Foam::endl;    
    
    //basisMesh.write();
    runTime.loop();
    runTime.write();

    
    // ---
    Info << nl << "To best visualise the results, load the mesh and extract all patches" << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

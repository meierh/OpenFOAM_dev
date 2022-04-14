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

#include "dynamicRefineFvMesh.H"
#include "cutCellFvMesh.H"
#include "Nurbs1D.H"
#include "UnitTestNurbs1D.H"
#include "UnitTestNurbs2D.H"
#include "KdTree.H"
#include "BsTree.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void UnitTests(int argc, char *argv[],Time& runTime)
{
    Info<<"---------------------------Unit Test---------------------------"<<endl;
    //TESTNURBS1D TestClass1;
    TESTNURBS2D TestClass2;

    //UnitTest_KdTree();
    //UnitTest_BsTree();
    //UnitTest_cutCellFvMesh(argc,argv,runTime);
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
    
    std::unique_ptr<scalarList> knots(new scalarList(12));
    (*knots)[0] = 0;    (*knots)[1] = 0;    (*knots)[2] = 0;    (*knots)[3] = 1;
    (*knots)[4] = 1;    (*knots)[5] = 2;    (*knots)[6] = 2;    (*knots)[7] = 3;
    (*knots)[8] = 3;    (*knots)[9] = 4;    (*knots)[10] = 4;   (*knots)[11] = 4;
    Info<<"Knoten"<<endl;
    
    //int testdegree = 2;
    
    std::unique_ptr<scalarList> weights(new scalarList(9));
    (*weights)[0] = 1;    (*weights)[1] = sqrt(2)/2;    (*weights)[2] = 1;
    (*weights)[3] = sqrt(2)/2;    (*weights)[4] = 1;    (*weights)[5] = sqrt(2)/2;
    (*weights)[6] = 1;    (*weights)[7] = sqrt(2)/2;    (*weights)[8] = 1;
    Info<<"Gewichte"<<endl;
    
    std::unique_ptr<List<Foam::vector>> controlPoints(new List<Foam::vector>(9));
    (*controlPoints)[0] = Foam::vector(1,0,0);    (*controlPoints)[1] = Foam::vector(1,1,0);    (*controlPoints)[2] = Foam::vector(0,1,0);
    (*controlPoints)[3] = Foam::vector(-1,1,0);   (*controlPoints)[4] = Foam::vector(-1,0,0);   (*controlPoints)[5] = Foam::vector(-1,-1,0);
    (*controlPoints)[6] = Foam::vector(0,-1,0);   (*controlPoints)[7] = Foam::vector(1,-1,0);   (*controlPoints)[8] = Foam::vector(1,0,0);
    Info<<"Kontrollpunkte"<<endl;
    
    runTime.loop();
    runTime.write();

    
    // ---
    Info << nl << "To best visualise the results, load the mesh and extract all patches" << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

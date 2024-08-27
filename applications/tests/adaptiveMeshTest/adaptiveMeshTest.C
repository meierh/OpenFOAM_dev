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
    Transient solver for incompressible, laminar flow of Newtonian fluids.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "pisoControl.H"
#include "fvModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"    
    Foam::pisoControl piso(mesh);
    #include "createFields.H"
    Foam::scalar cumulativeContErr = 0;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    Foam::Info<< "\nStarting time loop\n" << Foam::endl;
    Foam::label timeStep = 0;
    while (runTime.loop())
    {
        Foam::Info<< "Time = " << runTime.name() << Foam::nl << Foam::endl;

        Foam::Info<<"Topochanging:"<<mesh.topoChanging()<<Foam::endl;

        for(Foam::label i=0; i<refine.size(); i++)
        {
            refine[i] = 1;
        }
        Foam::Info<<"Refine field size:"<<refine.size()<<Foam::endl;
        bool refine1 = mesh.update();
        Foam::Info<<"refine 1:"<<refine1<<Foam::endl;
        Foam::Info<<"Topochanging:"<<mesh.topoChanging()<<Foam::endl;
        runTime.write();


        for(Foam::label i=0; i<refine.size(); i++)
        {
            refine[i] = 0;
        }
        Foam::Info<<"Refine field size:"<<refine.size()<<Foam::endl;
        bool refine2 = mesh.update();
        Foam::Info<<"refine 2:"<<refine2<<Foam::endl;
        Foam::Info<<"Topochanging:"<<mesh.topoChanging()<<Foam::endl;
        runTime.write();

        for(Foam::label i=0; i<refine.size(); i++)
        {
            refine[i] = -1;
        }
        Foam::Info<<"Refine field size:"<<refine.size()<<Foam::endl;
        bool refine3 = mesh.update();
        Foam::Info<<"refine 3:"<<refine3<<Foam::endl;
        Foam::Info<<"Topochanging:"<<mesh.topoChanging()<<Foam::endl;
        runTime.write();

        timeStep++;
        runTime.write();

        Foam::Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << Foam::nl << Foam::endl;
    }
    Foam::Info<< "End\n" << Foam::endl;

    return 0;
}


// ************************************************************************* //

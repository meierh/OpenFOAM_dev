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

#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "adjustPhi.H"
#include "fvcMeshPhi.H"
#include "fvcFlux.H"
#include "fvcDdt.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvmLaplacian.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"    
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    Foam::Info<< "\nStarting time loop\n" << Foam::endl;
    Foam::label timeStep = 0;
    while (runTime.loop())
    {
        Foam::Info<< "Time = " << runTime.name() << Foam::nl << Foam::endl;

        Foam::Info<<"Topochanging:"<<mesh.topoChanging()<<Foam::endl;

        Foam::dimensionSet dim(0,2,-2,0,0,0,0);
        Foam::dimensioned<Foam::scalar> val("test",dim,0);

        val.value() = 1;
        testField = val;
        for(Foam::label i=0; i<refine.size(); i++)
        {
            refine[i] = 1;
        }
        Foam::Info<<"refine:"<<refine.primitiveField()<<Foam::endl;
        Foam::Info<<"Refine field size:"<<refine.size()<<Foam::endl;
        bool refine1 = mesh.update();
        Foam::Info<<"refine 1:"<<refine1<<Foam::endl;
        Foam::Info<<"Refine field size:"<<refine.size()<<Foam::endl;        Foam::Info<<"Topochanging:"<<mesh.topoChanging()<<Foam::endl<<Foam::endl;
        runTime.write();

        val.value() = 2;
        testField = val;
        for(Foam::label i=0; i<refine.size(); i++)
        {
            refine[i] = 0;
        }
        Foam::Info<<"refine:"<<refine.primitiveField()<<Foam::endl;
        Foam::Info<<"Refine field size:"<<refine.size()<<Foam::endl;
        bool refine2 = mesh.update();
        Foam::Info<<"Refine field size:"<<refine.size()<<Foam::endl;
        Foam::Info<<"refine 2:"<<refine2<<Foam::endl;
        Foam::Info<<"Topochanging:"<<mesh.topoChanging()<<Foam::endl<<Foam::endl;
        runTime.write();

        val.value() = 3;
        testField = val;
        for(Foam::label i=0; i<refine.size(); i++)
        {
            refine[i] = -1;
        }
        Foam::Info<<"field:"<<refine.primitiveField()<<Foam::endl;
        Foam::Info<<"Refine field size:"<<refine.size()<<Foam::endl;
        bool refine3 = mesh.update();
        Foam::Info<<"Refine field size:"<<refine.size()<<Foam::endl;
        Foam::Info<<"refine 3:"<<refine3<<Foam::endl;
        Foam::Info<<"Topochanging:"<<mesh.topoChanging()<<Foam::endl<<Foam::endl;
        runTime.write();

        timeStep++;
        runTime.write();

        Foam::Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << Foam::nl << Foam::endl;
            
        Foam::Info<<"mesh.toc():"<<mesh.toc()<<Foam::endl;
    }
    const Foam::fvMeshTopoChanger& topoCh = mesh.topoChanger();
    const Foam::fvMeshTopoChanger* topoChPtr = &topoCh;
        
    Foam::Info<< "End\n" << Foam::endl;

    return 0;
}


// ************************************************************************* //

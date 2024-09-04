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
 #include "pisoControl.H"
 #include "pressureReference.H"
 #include "findRefCell.H"
 #include "constrainPressure.H"
 #include "constrainHbyA.H"
 #include "adjustPhi.H"

 #include "fvcDdt.H"
 #include "fvcGrad.H"
 #include "fvcFlux.H"

 #include "fvmDdt.H"
 #include "fvmDiv.H"
 #include "fvmLaplacian.H"

#include "CrossSectionStructure.H"
#include "StaticVelocityPressureAction.H"
#include <memory>

using namespace Foam;
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Info<<" ------ icoFoamStaticRod_substrU ------ "<<Foam::endl;

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"    
    pisoControl piso(mesh);
    #include "createFields.H"
    #include "initContinuityErrs.H"


    IOobject testIO("pxx","structure",runTime,IOobject::MUST_READ,IOobject::NO_WRITE);
    Info<<"testIO:"<<testIO.filePath("",true).empty()<<Foam::endl;
    IOdictionary dict(testIO);
    Info<<dict.name()<<Foam::endl;

    /*
    IOobject obj
    (
        "p",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );
    Info<<"fileName:"<<obj.filePath("",true)<<Foam::endl;

    volScalarField p_
    (
        IOobject
        (
            "pxx",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    */

    ITstream refineInterval = dict.lookup("refineInterval");
    Info<<refineInterval.name()<<Foam::endl;
    Info<<refineInterval<<Foam::endl;
    token refineIntervalToken;
    refineInterval.read(refineIntervalToken);
    scalar num = refineIntervalToken.scalarToken();
    Info<<"num:"<<num<<Foam::endl;

    ITstream rodFile = dict.lookup("rodType");
    Info<<rodFile.name()<<Foam::endl;
    Info<<rodFile<<Foam::endl;
    token rodType;
    rodFile.read(rodType);
    std::string rodTypeStr = rodType.stringToken();
    Info<<"rodTypeStr:"<<rodTypeStr<<Foam::endl;

    Info<<"rodType:"<<rodType<<Foam::endl;

    const dictionary& crossSecDict = dict.subDict("boundary");
    Info<<crossSecDict.tokens()<<Foam::endl;
    List<keyType> keys = crossSecDict.keys();
    for(keyType tok : keys)
    {
        const dictionary& boundaryDict = crossSecDict.subDict(tok);
        ITstream rodFile = boundaryDict.lookup("faces");
        while(rodFile.nRemainingTokens()>0)
        {
            token tok;
            rodFile.read(tok);
            Info<<tok<<Foam::endl;
        }
        break;
    }



    FatalErrorInFunction<<"Temp stop"<<exit(FatalError);

    CrossSection circle(0.05);
    CrossSectionStructure structure(mesh,{circle});

    StaticVelocityPressureAction U_Interaction(mesh,structure,U,Uf);

    std::string fileName = "rodForce.data";
    std::ofstream rodForceData;
    if(Pstream::master())
        rodForceData = std::ofstream(fileName);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    Info<< "\nStarting time loop\n" << endl;
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.name() << nl << endl;
        #include "CourantNo.H"

        bool converged = false;
        do
        {
            SolverPerformance<vector> solverResU;
            SolverPerformance<scalar> solverResP;

            fvVectorMatrix UEqn
            (
                fvm::ddt(U)
                + fvm::div(phi, U)
                - fvm::laplacian(nu, U)
            );
            if (piso.momentumPredictor())
            {
                solverResU = solve(UEqn == -fvc::grad(p) + Uf);
            }
            U_Interaction.interpolateFluidVelocityToMarkers();
            U_Interaction.computeCouplingForceOnMarkers();
            U_Interaction.computeRodForceMoment();
            U_Interaction.interpolateFluidForceField();

            U = U + runTime.deltaT()*Uf;

            while (piso.correct())
            {
                volScalarField rAU(1.0/UEqn.A());
                volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
                surfaceScalarField phiHbyA
                (
                    "phiHbyA",
                    fvc::flux(HbyA)
                    + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
                );
                adjustPhi(phiHbyA, U, p);

                // Update the pressure BCs to ensure flux consistency
                constrainPressure(p, U, phiHbyA, rAU);

                // Non-orthogonal pressure corrector loop
                while (piso.correctNonOrthogonal())
                {
                    // Pressure corrector
                    fvScalarMatrix pEqn
                    (
                        fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                    );
                    pEqn.setReference(pRefCell, pRefValue);
                    solverResP = pEqn.solve();

                    if (piso.finalNonOrthogonalIter())
                    {
                        phi = phiHbyA - pEqn.flux();
                    }
                }

                #include "continuityErrs.H"
                U = HbyA - rAU*fvc::grad(p);
                U.correctBoundaryConditions();
            }

            vector nIterationsU = solverResU.nIterations();
            converged = nIterationsU[0]>0?false:true;
            converged = nIterationsU[1]>0?false:true;
            converged = nIterationsU[2]>0?false:true;
            label nIterationsP = solverResP.nIterations();
            converged = nIterationsP>0?false:true;
        }
        while (!converged);

        FatalErrorInFunction<<"Temp stop"<<exit(FatalError);

        vector forces = U_Interaction.sumForces();
        Info<<"sum Force:"<<forces<<Foam::endl;
        if(Pstream::master())
        {
            rodForceData<<forces[0]<<" "<<forces[1]<<" "<<forces[2]<<std::endl;
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

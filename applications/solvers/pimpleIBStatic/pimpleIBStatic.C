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
    Info<<" ------ icoFoamStaticRod_addForce ------ "<<Foam::endl;

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"    
    pisoControl piso(mesh);
    #include "createFields.H"
    #include "initContinuityErrs.H"
    
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

        label iterationN = 0;
        bool converged = true;
        do
        {
            Info<<"-----------Iteration-------------"<<iterationN++<<Foam::endl;
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

            UEqn = fvVectorMatrix
            (
                fvm::ddt(U)
                + fvm::div(phi, U)
                - fvm::laplacian(nu, U)
            );
            if (piso.momentumPredictor())
            {
                solverResU = solve(UEqn == -fvc::grad(p) + Uf);
            }

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

                volVectorField UfbyA(Uf/UEqn.A());

                // Update the pressure BCs to ensure flux consistency
                constrainPressure(p, U, phiHbyA, rAU);

                // Non-orthogonal pressure corrector loop
                while (piso.correctNonOrthogonal())
                {
                    // Pressure corrector
                    fvScalarMatrix pEqn
                    (
                        fvm::laplacian(rAU, p) == fvc::div(phiHbyA) + fvc::div(UfbyA)
                    );
                    pEqn.setReference(pRefCell, pRefValue);
                    solverResP = pEqn.solve();

                    if (piso.finalNonOrthogonalIter())
                    {
                        phi = phiHbyA - pEqn.flux();
                    }
                }

                #include "continuityErrs.H"
                U = HbyA - rAU*fvc::grad(p) + UfbyA;
                U.correctBoundaryConditions();
            }

            converged = true;

            vector nIterationsU = solverResU.nIterations();
            converged &= nIterationsU[0]>0?false:true;
            converged &= nIterationsU[1]>0?false:true;
            converged &= nIterationsU[2]>0?false:true;
            label nIterationsP = solverResP.nIterations();
            Info<<"Iterations:"<<nIterationsU<<" / "<<nIterationsP<<Foam::endl;
            converged &= nIterationsP>0?false:true;

            if(iterationN>100)
                converged = true;
        }
        while (!converged);

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

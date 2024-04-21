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

#include "LineStructure.H"
#include "CrossSectionStructure.H"
#include "MarkerImplementation.H"
#include <memory>
#include "fvCFD.H"
#include "dynamicRefineFvMesh.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"    
    pisoControl piso(mesh);
    #include "createFields.H"
    #include "initContinuityErrs.H"

    /*
    scalar pi = constant::mathematical::pi;
    CrossSection crossSecN(0.1);   
    CrossSectionStructure structure(mesh,alpha,T,p,U,nu,{crossSecN});
    VelocityPressureForceInteraction Uf_Interaction(mesh,U,p,Uf);
    structure.transferMarkers(Uf_Interaction);
    */
    
    List<scalar> crossSecArea = {1};
    LineStructure structure(mesh,alpha,T,p,U,nu,crossSecArea);
    VelocityPressureForceInteraction Uf_Interaction(mesh,U,p,Uf);
    structure.transferMarkers(Uf_Interaction);
    
    Uf_Interaction.interpolateFluidVelocityToMarkers();
    Uf_Interaction.computeCouplingForceOnMarkers();
    Uf_Interaction.interpolateFluidForceField();
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    FatalErrorInFunction<<"Temp Stop"<<exit(FatalError);
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        // Momentum predictor
        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
        );
        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));
        }
        
        
        structure.transferMarkers(Uf_Interaction);
        Uf_Interaction.interpolateFluidVelocityToMarkers();
        Uf_Interaction.computeCouplingForceOnMarkers();
        Uf_Interaction.interpolateFluidForceField();
        //structure.transferMarkers(Uf_Interaction);
        
        /*
        U = U + Uf;
        */

        // --- PISO loop
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
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)// + fvm::div(Uf);
                );

                pEqn.setReference(pRefCell, pRefValue);

                pEqn.solve();

                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
        }

        solve
        (
            fvm::ddt(T)
            + fvm::div(phi,T)
            - fvm::laplacian(alpha,T)
        );
        
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

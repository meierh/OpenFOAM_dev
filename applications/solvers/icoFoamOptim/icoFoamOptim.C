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

#include <nlopt.hpp>
#include "LineStructure.H"
//#include "CrossSectionStructure.H"
#include <memory>
#include "fvCFD.H"
#include "dynamicRefineFvMesh.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

double objFunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
    if (!grad.empty()) {
        grad[0] = 0.0;
        grad[1] = 0.5 / std::sqrt(x[1]);
    }
    return std::sqrt(x[1]);
}

typedef struct {
    double a, b;
} my_constraint_data;

double constraintFunc(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
    my_constraint_data *d = reinterpret_cast<my_constraint_data*>(data);
    double a = d->a, b = d->b;
    if (!grad.empty()) {
        grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
        grad[1] = -1.0;
    }
    return ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
}

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"    
    pisoControl piso(mesh);
    #include "createFields.H"
    #include "initContinuityErrs.H"

    nlopt::opt opt(nlopt::LD_MMA, 2);
    std::vector<double> lb(2);
    lb[0] = -HUGE_VAL; lb[1] = 0;
    opt.set_lower_bounds(lb);
    opt.set_min_objective(objFunc, NULL);
    my_constraint_data data[2] = { {2,0}, {-1,1} };
    opt.add_inequality_constraint(constraintFunc, &data[0], 1e-8);
    opt.add_inequality_constraint(constraintFunc, &data[1], 1e-8);
    opt.set_xtol_rel(1e-4);
    std::vector<double> x(2);
    x[0] = 1.234; x[1] = 5.678;
    double minf;

    try
    {
        nlopt::result result = opt.optimize(x, minf);
        std::cout << "found minimum at f(" << x[0] << "," << x[1] << ") = "
            << std::setprecision(10) << minf << std::endl;
    }
    catch(std::exception &e)
    {
        std::cout << "nlopt failed: " << e.what() << std::endl;
    } 

    
    /*
    
    LineStructure structure(mesh,{0.1});
    
    VelocityPressureForceInteraction Uf_Interaction(mesh,structure,U,Uf);
    structure.connect(Uf_Interaction);
    
    FatalErrorInFunction<<"Temp Stop"<< exit(FatalError); 
    */
    
    /*
    VelocityPressureForceInteraction Uf_Interaction(mesh,U,Uf);
    structure.transferMarkers(Uf_Interaction);


    
    scalar pi = constant::mathematical::pi;
    CrossSection crossSecN(0.05);
    CrossSectionStructure structure(mesh,{crossSecN});

    VelocityPressureForceInteraction Uf_Interaction(mesh,U,Uf);
    structure.transferMarkers(Uf_Interaction);
    */

    /*
    scalar pi = constant::mathematical::pi;
    CrossSection crossSecN(0.1);   
    CrossSectionStructure structure(mesh,{crossSecN});
    VelocityPressureForceInteraction Uf_Interaction(mesh,U,Uf);
    structure.transferMarkers(Uf_Interaction);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
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
        
        Uf_Interaction.interpolateFluidVelocityToMarkers();
        Uf_Interaction.computeCouplingForceOnMarkers();
        Uf_Interaction.computeRodForceMoment();
        Uf_Interaction.interpolateFluidForceField();
        Info<<"sum Force:"<<Uf_Interaction.sumForces()<<Foam::endl;
        
        U = U + runTime.deltaT()*Uf;

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
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
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
        
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
    Info<< "End\n" << endl;
    */

    return 0;
}


// ************************************************************************* //

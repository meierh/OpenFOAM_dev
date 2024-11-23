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
    foamRun

Description
    Loads and executes an OpenFOAM solver module either specified by the
    optional \c solver entry in the \c controlDict or as a command-line
    argument.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient and steady simulations.

Usage
    \b foamRun [OPTION]

      - \par -solver <name>
        Solver name

      - \par -libs '(\"lib1.so\" ... \"libN.so\")'
        Specify the additional libraries loaded

    Example usage:
      - To run a \c rhoPimpleFoam case by specifying the solver on the
        command line:
        \verbatim
            foamRun -solver fluid
        \endverbatim

      - To update and run a \c rhoPimpleFoam case add the following entries to
        the controlDict:
        \verbatim
            application     foamRun;

            solver          fluid;
        \endverbatim
        then execute \c foamRun

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Optimizer.H"
#include "icoAdjointImmersedBoundary.H"
#include "pimpleSingleRegionControl.H"
#include "setDeltaT.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::solvers::icoAdjointImmersedBoundary::objectiveFunction createTotalPressureLoss()
{
    Foam::solvers::icoAdjointImmersedBoundary::objectiveFunction obj;
    
    // J = int_in (p+0.5*u²) dS - int_out (p+0.5*u²) dS
    
    /* 
     * dJdp = int_in 1 dS - int_out 1 dS
     * dJdu = int_in u dS - int_out u dS
     * dJdT = 0
     */
    
    obj.dJdp_InletWall = [](const icoAdjointVelocityInletWallBC& bc)
    {
        return Field<scalar>(bc.patch().size(),1);
    };
    obj.dJdu_uOutlet = [](const icoAdjointVelocityOutletBC& bc)
    {
        const fvPatchField<vector>& u = bc.patch().lookupPatchField<volVectorField,vector>("U");
        return Field<vector>(-u);
    };
    obj.dJdu_pOutlet = [](const icoAdjointPressureOutletBC& bc)
    {
        const fvPatchField<vector>& u = bc.patch().lookupPatchField<volVectorField,vector>("U");
        return Field<vector>(-u);
    };
    obj.J = [](Foam::solvers::icoAdjointImmersedBoundary const& domain)
    {
        scalar J=0;
        const fvBoundaryMesh& domainBCs = domain.mesh.boundary();
        const volScalarField& p = domain.p;
        const volScalarField::Boundary& p_boundary = p.boundaryField();
        const volVectorField& u = domain.U;
        const volVectorField::Boundary& u_boundary = u.boundaryField();
        
        const label inletPatchInd = domainBCs.findIndex("inlet");
        const fvPatch& inletPatch = domainBCs[inletPatchInd];
        const scalarField& inFaceMag = inletPatch.magSf();
        tmp<vectorField> inNormals = -1*inletPatch.nf();
        const fvPatchField<scalar>& pInlet = p_boundary[inletPatchInd];
        const fvPatchField<vector>& UInlet = u_boundary[inletPatchInd];
        scalarField u_minN_inlet = UInlet & inNormals.ref();
        scalarField abs_u = 0.5*(u_minN_inlet * u_minN_inlet);
        scalarField p_plus_abs_u = abs_u+pInlet;
        scalarField int_p_plus_abs_u = p_plus_abs_u*inFaceMag;
        for(scalar val : int_p_plus_abs_u)
            J += val;

        const label outletPatchInd = domainBCs.findIndex("outlet");
        const fvPatch& outletPatch = domainBCs[outletPatchInd];
        const scalarField& outFaceMag = outletPatch.magSf();
        tmp<vectorField> outNormals = outletPatch.nf();
        const fvPatchField<scalar>& pOutlet = p_boundary[outletPatchInd];
        const fvPatchField<vector>& UOutlet = u_boundary[outletPatchInd];
        scalarField u_minN_outlet = UOutlet & outNormals.ref();
        abs_u = 0.5*(u_minN_outlet * u_minN_outlet);
        p_plus_abs_u = abs_u+pOutlet;
        int_p_plus_abs_u = p_plus_abs_u*outFaceMag;
        for(scalar val : int_p_plus_abs_u)
            J -= val;
        
        return J;
    };   
    
    return obj;
}

double objFunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
    Foam::solvers::icoAdjointImmersedBoundary* solver = static_cast<Foam::solvers::icoAdjointImmersedBoundary*>(my_func_data);
    if(solver==nullptr)
        FatalErrorInFunction<<"Failure to cast to solver"<<exit(FatalError);
    
    if (!grad.empty()) {
        grad[0] = 2*x[0];
        grad[1] = 2*x[1];
    }
    return x[0]*x[0]+x[1]*x[1];
}

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // Instantiate the solver
    Foam::solvers::icoAdjointImmersedBoundary solver(mesh,runTime,createTotalPressureLoss());
    solver.SolveSteadyAdjoint();
  
    //label optimCoeffsNbr=1;
  
    //Optimizer
    //Optimizer optimSolver(objFunc,optimCoeffsNbr,&solver);

    // Create the outer PIMPLE loop and control structure
    //pimpleSingleRegionControl pimple(solver.pimple);

    // Set the initial time-step
    //setDeltaT(runTime, solver);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    /*
    Info<< nl << "Starting time loop\n" << endl;

    while (pimple.run(runTime))
    {
        // Update PIMPLE outer-loop parameters if changed
        pimple.read();

        solver.preSolve();

        // Adjust the time-step according to the solver maxDeltaT
        adjustDeltaT(runTime, solver);

        runTime++;

        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        // PIMPLE corrector loop
        while (pimple.loop())
        {
            solver.moveMesh();
            solver.motionCorrector();
            solver.fvModels().correct();
            solver.prePredictor();
            solver.momentumPredictor();
            solver.thermophysicalPredictor();
            solver.pressureCorrector();
            solver.postCorrector();
        }

        solver.postSolve();

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

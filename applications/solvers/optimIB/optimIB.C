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

objectiveFunction createTotalPressureLoss()
{
    objectiveFunction obj;
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

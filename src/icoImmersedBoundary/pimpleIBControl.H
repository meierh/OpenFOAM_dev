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

Class
    Foam::solvers::icoImmersedBoundary

Description
    Solver module for steady or transient turbulent flow of incompressible
    isothermal fluids with optional mesh motion and change.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient and steady simulations.

    Optional fvModels and fvConstraints are provided to enhance the simulation
    in many ways including adding various sources, constraining or limiting
    the solution.

    Reference:
    \verbatim
        Greenshields, C. J., & Weller, H. G. (2022).
        Notes on Computational Fluid Dynamics: General Principles.
        CFD Direct Ltd.: Reading, UK.
    \endverbatim

SourceFiles
    icoImmersedBoundary.C

See also
    Foam::solvers::fluidSolver
    Foam::solvers::isothermalFluid

\*---------------------------------------------------------------------------*/

#ifndef pimpleIBControl_H
#define pimpleIBControl_H

#include "incompressibleFluid.H"
#include "CrossSectionStructure.H"
#include "StaticVelocityPressureAction.H"
#include "ForcedMovementVelocityPressureAction.H"
#include "FixedTemperatureAction.H"
#include "MeshRefiner.H"

#include "pimpleSingleRegionControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class incompressibleFluid Declaration
\*---------------------------------------------------------------------------*/

class pimpleIBControl : public pimpleSingleRegionControl
{
    public:
        pimpleIBControl(pimpleNoLoopControl& pimple,Time& runTime);

        bool loop()=delete;
        
        bool runNonStop(Time& time);
        bool run(Time& time);
        
        bool momentumLoopOrig();
        bool momentumLoop();
        bool momentumPredictor();
        bool consistent();
        label nCorrPiso();
        bool correct();
        bool correctNonOrthogonal();
        bool finalNonOrthogonalIter();

        bool temperatureLoop();
        
        void setVelocityPerformance(SolverPerformance<vector>* uEqn);
        void setPressurePerformance(SolverPerformance<scalar>* pEqn);
        void setTemperaturePerformance(SolverPerformance<scalar>* tEqn);
        
        bool isSteady(){return steady;}

        void reset();

    protected:
        pimpleNoLoopControl& pimple;        
        void readFromFvSolution();
        Time& runTime;
        bool steady;
        label timeIteration;
        
    private:
        // Momentum system
        SolverPerformance<vector>* velocityEqns=nullptr;
        scalar velocityTolerance;
        SolverPerformance<scalar>* pressureEqns=nullptr;
        scalar pressureTolerance;
        label momentumInnerIteration;
        label minMomentumIterations=1;
        scalar momentumIterTolFrac=1;
        
        // Temperature system
        SolverPerformance<scalar>* temperatureEqns;
        scalar temperatureTolerance;
        bool temperatureToleranceSet = false;
        label temperatureInnerIteration;
        label minTemperatureIterations=1;
        scalar temperatureIterTolFrac=1;
        bool temperatureUsed = false;
        
        //Full system
        label minOuterIterations=1; // nOuterCorrectors
        bool advancedInitialConvergence = true;
        scalar convergenceDelay;
        bool convergenceDelaySet;

        bool resetActive = false;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

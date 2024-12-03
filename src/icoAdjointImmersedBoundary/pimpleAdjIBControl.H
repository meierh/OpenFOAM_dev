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

#ifndef pimpleAdjIBControl_H
#define pimpleAdjIBControl_H

#include "incompressibleFluid.H"
#include "CrossSectionStructure.H"
#include "StaticVelocityPressureAction.H"
#include "ForcedMovementVelocityPressureAction.H"
#include "FixedTemperatureAction.H"
#include "MeshRefiner.H"

#include "pimpleIBControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class incompressibleFluid Declaration
\*---------------------------------------------------------------------------*/

class pimpleAdjIBControl : public pimpleIBControl
{
    public:
        pimpleAdjIBControl(pimpleNoLoopControl& pimple,Time& runTime);
        pimpleAdjIBControl(const pimpleAdjIBControl&)=delete;
        //virtual ~pimpleIBControl() {~pimpleSingleRegionControl::pimpleSingleRegionControl();}

        bool loop()=delete;
        
        bool run(Time& time);
        
        bool adjMomentumLoop();
        bool adjMomentumPredictor();
        bool adjConsistent();
        label adjNCorrPiso();
        bool adjCorrect();
        bool adjCorrectNonOrthogonal();
        bool adjFinalNonOrthogonalIter();
        
        bool adjTemperatureLoop();
        
        void setAdjVelocityPerformance(SolverPerformance<vector>* adjUEqn);
        void setAdjPressurePerformance(SolverPerformance<scalar>* adjPEqn);
        void setAdjTemperaturePerformance(SolverPerformance<scalar>* adjTEqn);
        
    private:
        void readFromFvSolution();
        
        // Momentum system
        SolverPerformance<vector>* adjVelocityEqns=nullptr;
        scalar adjVelocityTolerance;
        SolverPerformance<scalar>* adjPressureEqns=nullptr;
        scalar adjPressureTolerance;
        scalar adjMomentumTolerance;
        label momentumIteration = 0;
        
        // Temperature system
        SolverPerformance<scalar>* adjTemperatureEqns;
        scalar adjTemperatureTolerance;
        bool adjTemperatureToleranceSet = false;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

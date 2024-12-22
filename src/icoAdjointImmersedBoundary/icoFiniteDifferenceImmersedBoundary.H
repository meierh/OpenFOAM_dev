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

#ifndef icoFiniteDifferenceImmersedBoundary_H
#define icoFiniteDifferenceImmersedBoundary_H

#include "SetupOpenFOAM.H"
#include "icoAdjointImmersedBoundary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{    
/*---------------------------------------------------------------------------*\
                          Class incompressibleFluid Declaration
\*---------------------------------------------------------------------------*/

class icoFiniteDifferenceImmersedBoundary
{
    public:
    // Constructors
        //- Construct from region mesh
        icoFiniteDifferenceImmersedBoundary
        (
            argList& args,
            Parameter para,
            std::vector<scalar> percFD,
            objectiveFunction obj=createEmptyObjective()
        );
        //- Disallow default bitwise assignment
        void operator=(const icoAdjointImmersedBoundary&) = delete;
    
    // High level functions
        void Solve();
        
    protected:
        argList& args;
        scalar parameterIniValue;
        Parameter para;
        std::vector<scalar> epsilons;
        std::function<scalar(const icoAdjointImmersedBoundary&)> J;
        //std::unique_ptr<icoAdjointImmersedBoundary> icoSolver;
        std::unique_ptr<std::ofstream> recordFDFile;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

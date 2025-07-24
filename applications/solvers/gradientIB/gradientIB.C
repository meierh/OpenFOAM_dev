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
#include "icoAdjointImmersedBoundary.H"

#include "pimpleSingleRegionControl.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Pout<<"gradientIB start"<<Foam::endl;
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Parameter radPara(CrossSectionCoeffReference(0,0,0));

    // Instantiate the solver
    Foam::solvers::icoAdjointImmersedBoundary solver(mesh,runTime,{radPara});
    solver.SolveSteadyAdjoint();

    if(Pstream::master())
    {
        const std::vector<std::pair<Parameter,scalar>>& parameterGradient = solver.getParameterGradient();
        std::ofstream recordGradFile("gradRecords");
        for(const std::pair<Parameter,scalar>& singleParameter : parameterGradient)
        {
            Pout<<singleParameter.first.to_string()<<": "<<singleParameter.second<<Foam::nl;
            recordGradFile<<singleParameter.first.to_string()<<": "<<singleParameter.second<<std::endl;
        }
    }

    Info<<"gradientIB done"<<Foam::nl;

    return 0;
}



// ************************************************************************* //

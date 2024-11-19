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

\*---------------------------------------------------------------------------*/

#include "setDeltaT.H"
#include "icoImmersedBoundary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::setDeltaT(Time& runTime, const solver& solver)
{
    if
    (
        runTime.timeIndex() == 0
     && runTime.controlDict().lookupOrDefault("adjustTimeStep", false)
     && solver.transient()
    )
    {
        const scalar deltaT =
            min(solver.maxDeltaT(), runTime.functionObjects().maxDeltaT());

        if (deltaT < rootVGreat)
        {
            runTime.setDeltaT(min(runTime.deltaTValue(), deltaT));
        }
    }
}


void Foam::adjustDeltaT(Time& runTime, const solver& solver)
{
    // Update the time-step limited by the solver maxDeltaT
    if
    (
        runTime.controlDict().lookupOrDefault("adjustTimeStep", false)
     && solver.transient()
    )
    {
        const scalar deltaTNew = min(solver.maxDeltaT(), runTime.functionObjects().maxDeltaT());
        const Foam::solvers::icoImmersedBoundary* solIcoIB = dynamic_cast<const Foam::solvers::icoImmersedBoundary*>(&solver);
        if(solIcoIB==nullptr)
            FatalErrorInFunction<<"Failed casting to icoImmersedBoundary"<<exit(FatalError);
        
        const scalar minDeltaTLim = solIcoIB->minDeltaT();
        if(solver.maxDeltaT()<minDeltaTLim)
            FatalErrorInFunction<<"maxDeltaT is smaller than minDeltaT"<<exit(FatalError);
        const scalar deltaT = (deltaTNew<minDeltaTLim) ? minDeltaTLim : deltaTNew;

        if (deltaT < rootVGreat)
        {
            runTime.setDeltaT
            (
                min(solver::deltaTFactor*runTime.deltaTValue(), deltaT)
            );
            Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
        }
    }
}


// ************************************************************************* //

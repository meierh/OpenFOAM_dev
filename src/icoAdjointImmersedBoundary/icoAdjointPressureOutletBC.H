/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::icoAdjointPressureOutletBC

Description


SourceFiles
    icoAdjointPressureOutletBC.C

\*---------------------------------------------------------------------------*/

#ifndef icoAdjointPressureOutletBC_H
#define icoAdjointPressureOutletBC_H

#include "fixedValueFvPatchField.H"
#include <functional>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class icoAdjointPressureOutletBC Declaration
\*---------------------------------------------------------------------------*/

class icoAdjointPressureOutletBC : public fixedValueFvPatchField<scalar>
{
    public:
        //- Runtime type information
        TypeName("icoAdjointPressureOutlet");

    // Constructors
        //- Construct from patch, internal field and dictionary
        icoAdjointPressureOutletBC
        (
            const fvPatch&,
            const DimensionedField<scalar,volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given icoAdjointPressureOutletBC
        //  onto a new patch
        icoAdjointPressureOutletBC
        (
            const icoAdjointPressureOutletBC&,
            const fvPatch&,
            const DimensionedField<scalar,volMesh>&,
            const fieldMapper&
        );

        //- Construct as copy setting internal field reference
        icoAdjointPressureOutletBC
        (
            const icoAdjointPressureOutletBC&,
            const DimensionedField<scalar,volMesh>&
        );

    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();
        
        
        void set_dJdu_Outlet(std::function<Field<vector>(const icoAdjointPressureOutletBC&)> expr);
        void set_nu(dimensionedScalar nu);
        void set_temperatureUsed(bool temperatureUsed);

    private:
        std::function<Field<vector>(const icoAdjointPressureOutletBC&)> dJdu_Outlet;
        
        bool nu_set=false;
        dimensionedScalar nu;
        
        bool temperatureUsed_set=false;
        bool temperatureUsed;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

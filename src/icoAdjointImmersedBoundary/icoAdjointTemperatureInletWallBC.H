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
    Foam::icoAdjointVelocityOutletBC

Description

SourceFiles
    icoAdjointVelocityOutletBC.C

\*---------------------------------------------------------------------------*/

#ifndef icoAdjointTemperatureInletWallBC_H
#define icoAdjointTemperatureInletWallBC_H

#include "fixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class adjointOutletVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class icoAdjointTemperatureInletWallBC : public fixedValueFvPatchField<scalar>
{
    public:
        //- Runtime type information
        TypeName("icoAdjointTemperatureInletWall");

        // Constructors
        //- Construct from patch, internal field and dictionary
        icoAdjointTemperatureInletWallBC
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given icoAdjointTemperatureInletWallBC
        //  onto a new patch
        icoAdjointTemperatureInletWallBC
        (
            const icoAdjointTemperatureInletWallBC&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Construct as copy setting internal field reference
        icoAdjointTemperatureInletWallBC
        (
            const icoAdjointTemperatureInletWallBC&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new icoAdjointTemperatureInletWallBC(*this, iF)
            );
        }

    // Member Functions

        //- Update the coefficients associated with the patch field
        //virtual void updateCoeffs();

        //- Write
        //virtual void write(Ostream&) const;
        
    private:
        const dictionary& setZeroValue(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

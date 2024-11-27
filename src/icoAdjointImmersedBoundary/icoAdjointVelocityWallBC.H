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

#ifndef icoAdjointVelocityWallBC_H
#define icoAdjointVelocityWallBC_H

#include "fixedValueFvPatchField.H"
#include <functional>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class adjointOutletVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class icoAdjointVelocityWallBC : public fixedValueFvPatchField<vector>
{
    public:
        //- Runtime type information
        TypeName("icoAdjointVelocityWall");

        // Constructors
        //- Construct from patch, internal field and dictionary
        icoAdjointVelocityWallBC
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given icoAdjointVelocityWallBC
        //  onto a new patch
        icoAdjointVelocityWallBC
        (
            const icoAdjointVelocityWallBC&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fieldMapper&
        );

        //- Construct as copy setting internal field reference
        icoAdjointVelocityWallBC
        (
            const icoAdjointVelocityWallBC&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new icoAdjointVelocityWallBC(*this, iF)
            );
        }

    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();
        
        void set_dJdp_Wall(std::function<Field<scalar>(const icoAdjointVelocityWallBC&)> expr);

    private:
        std::function<Field<scalar>(const icoAdjointVelocityWallBC&)> dJdp_Wall;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

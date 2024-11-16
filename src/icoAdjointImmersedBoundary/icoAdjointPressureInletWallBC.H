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
    Foam::icoAdjointPressureInletWallBC

Description
    a phi + b dphi/dn = c

SourceFiles
    icoAdjointPressureInletWallBC.C

\*---------------------------------------------------------------------------*/

#ifndef icoAdjointPressureInletWallBC_H
#define icoAdjointPressureInletWallBC_H

#include "uniformFixedGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class icoAdjointPressureInletWallBC Declaration
\*---------------------------------------------------------------------------*/

class icoAdjointPressureInletWallBC : public uniformFixedGradientFvPatchField<scalar>
{
    public:
        //- Runtime type information
        TypeName("icoAdjointPressureInletWall");

    // Constructors
        //- Construct from patch and internal field
        /*
        icoAdjointPressureInletWallBC
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const Field<scalar>&
        );
        */

        //- Construct from patch, internal field and dictionary
        icoAdjointPressureInletWallBC
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given icoAdjointPressureInletWallBC
        //  onto a new patch
        
        icoAdjointPressureInletWallBC
        (
            const icoAdjointPressureInletWallBC&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );
        

        //- Disallow copy without setting internal field reference
        icoAdjointPressureInletWallBC
        (
            const icoAdjointPressureInletWallBC&
        ) = delete;

        //- Copy constructor setting internal field reference
        icoAdjointPressureInletWallBC
        (
            const icoAdjointPressureInletWallBC&,
            const DimensionedField<scalar, volMesh>&
        );
        

        //- Construct and return a clone setting internal field reference
        /*
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar>>
            (
                new icoAdjointPressureInletWallBC(*this, iF)
            );
        }
        */
        

    // Member Functions

        //- Update the coefficients associated with the patch field
        //virtual void updateCoeffs();

        //- Write
        //virtual void write(Ostream&) const;       
        
    private:
        const dictionary& setZeroGrad(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

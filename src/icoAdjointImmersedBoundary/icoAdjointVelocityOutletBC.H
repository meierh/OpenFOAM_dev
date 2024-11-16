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

#ifndef icoAdjointVelocityOutletBC_H
#define icoAdjointVelocityOutletBC_H

#include "robinFvPatchField.H"
#include <functional>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class adjointOutletVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class icoAdjointVelocityOutletBC : public robinFvVectorPatchField
{
    public:

        //- Runtime type information
        TypeName("icoAdjointVelocityOutlet");

    // Constructors
        //- Construct from patch, internal field and dictionary
        icoAdjointVelocityOutletBC
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given icoAdjointVelocityOutletBC
        //  onto a new patch
        icoAdjointVelocityOutletBC
        (
            const icoAdjointVelocityOutletBC&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fieldMapper&
        );

        //- Construct as copy setting internal field reference
        icoAdjointVelocityOutletBC
        (
            const icoAdjointVelocityOutletBC&,
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
                new icoAdjointVelocityOutletBC(*this, iF)
            );
        }

    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();
        
        void set_dJdu_Outlet(std::function<Field<vector>(const icoAdjointVelocityOutletBC&)> expr);
        void set_nu(dimensionedScalar nu);

    private:
        std::function<Field<vector>(const icoAdjointVelocityOutletBC&)> dJdu_Outlet;
        
        bool nu_set=false;
        dimensionedScalar nu;
        
        void write_a() override;
        void write_b() override;
        void write_c() override;        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

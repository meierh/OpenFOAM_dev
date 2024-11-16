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

#ifndef icoAdjointTemperatureOutletBC_H
#define icoAdjointTemperatureOutletBC_H

#include "robinFvPatchField.H"
#include <functional>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class adjointOutletVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class icoAdjointTemperatureOutletBC : public robinFvScalarPatchField
{
    public:
        //- Runtime type information
        TypeName("icoAdjointTemperatureOutlet");

    // Constructors
        //- Construct from patch, internal field and dictionary
        icoAdjointTemperatureOutletBC
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given icoAdjointTemperatureOutletBC
        //  onto a new patch
        icoAdjointTemperatureOutletBC
        (
            const icoAdjointTemperatureOutletBC&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Construct as copy setting internal field reference
        icoAdjointTemperatureOutletBC
        (
            const icoAdjointTemperatureOutletBC&,
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
                new icoAdjointTemperatureOutletBC(*this, iF)
            );
        }

    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        void set_dJdT_Outlet(std::function<Field<scalar>(const icoAdjointTemperatureOutletBC&)> expr);
        void set_alpha(dimensionedScalar alpha);

    private:
        std::function<Field<scalar>(const icoAdjointTemperatureOutletBC&)> dJdT_Outlet;
        
        bool alpha_set=false;
        dimensionedScalar alpha;
        
        void write_a() override;
        void write_b() override;
        void write_c() override;        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

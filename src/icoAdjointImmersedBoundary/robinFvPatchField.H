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
    Foam::robinFvPatchField

Description
    a phi + b dphi/dn = c

SourceFiles
    robinFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef robinFvPatchField_H
#define robinFvPatchField_H

#include "fixedValueFvPatchField.H"

#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class robinFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class robinFvPatchField : public fixedValueFvPatchField<Type>
{
    public:

    // Constructors
        //- Construct from patch, internal field and dictionary
        robinFvPatchField
        (
            const fvPatch& p,
            const DimensionedField<Type, volMesh>& iF,
            const dictionary& dict
        ):
        fixedValueFvPatchField<Type>(p,iF,dict,false)
        {}

        //- Construct by mapping given robinFvPatchField
        //  onto a new patch
        robinFvPatchField
        (
            const robinFvPatchField& ptf,
            const fvPatch& p,
            const DimensionedField<Type, volMesh>& iF,
            const fieldMapper& mapper
        ):
        fixedValueFvPatchField<Type>(ptf,p,iF,mapper)
        {}

        //- Construct as copy setting internal field reference
        robinFvPatchField
        (
            const robinFvPatchField& pivpvf,
            const DimensionedField<Type, volMesh>& iF
        ):
        fixedValueFvPatchField<Type>(pivpvf,iF)
        {}
        
        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs()
        {
            if (updated())
            {
                return;
            }
            
            tmp<Field<Type>> cellValues = patch().patchInternalField();
            
            scalarField dx = computeFaceCellDist();
            
            scalarField b_div_dx = b/dx;
            
            scalarField a_min_b_div_dx = a-b_div_dx;
            scalarField inv_a_min_b_div_dx = 1/a_min_b_div_dx;
            
            Field<Type> b_div_dx_mul_val = b_div_dx*cellValues;
            Field<Type> c_min_b_div_dx_mul_val = c-b_div_dx_mul_val;
            
            Field<Type>::operator=(inv_a_min_b_div_dx * c_min_b_div_dx_mul_val);
            
            fixedValueFvPatchField<vector>::updateCoeffs(); // sets updated_ to true
        }
        
        scalarField computeFaceCellDist()
        {
            const vectorField& faceCentres = patch().Cf();
            tmp<vectorField> cellCentres = patch().Cn();
            tmp<vectorField> faceNormals = patch().nf();
            vectorField centreToFace = faceCentres-cellCentres;
            scalarField centreToFaceDist = centreToFace & faceNormals;
            return centreToFaceDist;
        }

        //- Write
        //virtual void write(Ostream&) const;
        
    protected:
        virtual void write_a()=0;
        virtual void write_b()=0;
        virtual void write_c()=0;
        
        Field<scalar> a;
        Field<scalar> b;
        Field<Type> c;
};

using robinFvScalarPatchField = robinFvPatchField<scalar>;
using robinFvVectorPatchField = robinFvPatchField<vector>;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "robinFvPatchField.T.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

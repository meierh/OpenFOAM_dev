#include "icoAdjointTemperatureInletWallBC.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::icoAdjointTemperatureInletWallBC::icoAdjointTemperatureInletWallBC
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:fixedValueFvPatchField<scalar>(p, iF, dict, false)
{}

Foam::icoAdjointTemperatureInletWallBC::icoAdjointTemperatureInletWallBC
(
    const icoAdjointTemperatureInletWallBC& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:fixedValueFvPatchField<scalar>(ptf, p, iF, mapper)
{}


Foam::icoAdjointTemperatureInletWallBC::icoAdjointTemperatureInletWallBC
(
    const icoAdjointTemperatureInletWallBC& pivpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:fixedValueFvPatchField<scalar>(pivpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
const Foam::dictionary& Foam::icoAdjointTemperatureInletWallBC::setZeroValue(const dictionary& dict)
{
    dictionary& nonContDict = const_cast<dictionary&>(dict);
    nonContDict.set("uniform",0);
    return dict;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        icoAdjointTemperatureInletWallBC
    );
}


// ************************************************************************* //

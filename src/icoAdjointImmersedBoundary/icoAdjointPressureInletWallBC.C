#include "icoAdjointPressureInletWallBC.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
/*
Foam::icoAdjointPressureInletWallBC::icoAdjointPressureInletWallBC
(
    const fvPatch& p,
    const DimensionedField<scalar,volMesh>& iF,
    const Field<scalar>& fld
)
:uniformFixedGradientFvPatchField<scalar>(p, iF, fld)
{}
*/

Foam::icoAdjointPressureInletWallBC::icoAdjointPressureInletWallBC
(
    const fvPatch& p,
    const DimensionedField<scalar,volMesh>& iF,
    const dictionary& dict
)
:uniformFixedGradientFvPatchField<scalar>(p, iF, setZeroGrad(dict))
{}

Foam::icoAdjointPressureInletWallBC::icoAdjointPressureInletWallBC
(
    const icoAdjointPressureInletWallBC& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:uniformFixedGradientFvPatchField<scalar>(ptf, p, iF, mapper)
{}

Foam::icoAdjointPressureInletWallBC::icoAdjointPressureInletWallBC
(
    const icoAdjointPressureInletWallBC& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:uniformFixedGradientFvPatchField<scalar>(ptf, iF)
{}

const Foam::dictionary& Foam::icoAdjointPressureInletWallBC::setZeroGrad(const dictionary& dict)
{
    dictionary& nonContDict = const_cast<dictionary&>(dict);
    nonContDict.set("uniformGradient",0);
    return dict;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        icoAdjointPressureInletWallBC
    );
}

// ************************************************************************* //

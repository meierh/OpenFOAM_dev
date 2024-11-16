#include "icoAdjointPressureOutletBC.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::icoAdjointPressureOutletBC::icoAdjointPressureOutletBC
(
    const fvPatch& p,
    const DimensionedField<scalar,volMesh>& iF,
    const dictionary& dict
)
:fixedValueFvPatchField<scalar>(p, iF, dict, false)
{}

Foam::icoAdjointPressureOutletBC::icoAdjointPressureOutletBC
(
    const icoAdjointPressureOutletBC& ptf,
    const fvPatch& p,
    const DimensionedField<scalar,volMesh>& iF,
    const fieldMapper& mapper
)
:fixedValueFvPatchField<scalar>(ptf, p, iF, mapper)
{}

Foam::icoAdjointPressureOutletBC::icoAdjointPressureOutletBC
(
    const icoAdjointPressureOutletBC& pivpvf,
    const DimensionedField<scalar,volMesh>& iF
)
:fixedValueFvPatchField<scalar>(pivpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void Foam::icoAdjointPressureOutletBC::updateCoeffs()
{
    if(fixedValueFvPatchField<scalar>::updated())
    {
        return;
    }

}

void Foam::icoAdjointPressureOutletBC::set_dJdu_Outlet
(
    std::function<Field<vector>(const icoAdjointPressureOutletBC&)> expr
)
{
    Info<<"Set dJdu_pOutlet"<<Foam::nl;
    dJdu_Outlet = expr;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        icoAdjointPressureOutletBC
    );
}

// ************************************************************************* //

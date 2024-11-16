#include "icoAdjointVelocityInletWallBC.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::icoAdjointVelocityInletWallBC::icoAdjointVelocityInletWallBC
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:fixedValueFvPatchField<vector>(p, iF, dict, false)
{}

Foam::icoAdjointVelocityInletWallBC::icoAdjointVelocityInletWallBC
(
    const icoAdjointVelocityInletWallBC& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fieldMapper& mapper
)
:fixedValueFvPatchField<vector>(ptf, p, iF, mapper)
{}

Foam::icoAdjointVelocityInletWallBC::icoAdjointVelocityInletWallBC
(
    const icoAdjointVelocityInletWallBC& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:fixedValueFvPatchField<vector>(pivpvf, iF)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void Foam::icoAdjointVelocityInletWallBC::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    if(!dJdp_InletWall)
        FatalErrorInFunction<<"dJdp_InletWall not set"<<exit(FatalError);    
    scalarField dJdp = dJdp_InletWall(*this);
    vectorField n = patch().nf();
    vectorField::operator=(n*-1*dJdp);
    
    fixedValueFvPatchField<vector>::updateCoeffs(); // sets updated_ to true
}

void Foam::icoAdjointVelocityInletWallBC::set_dJdp_InletWall
(
    std::function<Field<scalar>(const icoAdjointVelocityInletWallBC&)> expr
)
{
    Info<<"Set dJdp_InletWall"<<Foam::nl;
    dJdp_InletWall = expr;
}

/*
void Foam::icoAdjointVelocityInletWallBC::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "value", *this);
}
*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        icoAdjointVelocityInletWallBC
    );
}


// ************************************************************************* //

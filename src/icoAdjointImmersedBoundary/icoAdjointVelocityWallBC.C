#include "icoAdjointVelocityWallBC.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::icoAdjointVelocityWallBC::icoAdjointVelocityWallBC
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:fixedValueFvPatchField<vector>(p, iF, dict, false)
{}

Foam::icoAdjointVelocityWallBC::icoAdjointVelocityWallBC
(
    const icoAdjointVelocityWallBC& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fieldMapper& mapper
)
:fixedValueFvPatchField<vector>(ptf, p, iF, mapper)
{}

Foam::icoAdjointVelocityWallBC::icoAdjointVelocityWallBC
(
    const icoAdjointVelocityWallBC& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:fixedValueFvPatchField<vector>(pivpvf, iF)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void Foam::icoAdjointVelocityWallBC::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    if(!dJdp_Wall)
        FatalErrorInFunction<<"dJdp_Wall not set"<<exit(FatalError);    
    scalarField dJdp = dJdp_Wall(*this);
    vectorField n = patch().nf();
    vectorField::operator=(n*-1*dJdp);
    
    fixedValueFvPatchField<vector>::updateCoeffs(); // sets updated_ to true

    Info<<"---------------------------------------------"<<Foam::nl;
    Info<<"| icoAdjointVelocityWallBC::updateCoeffs done"<<Foam::nl;
    vector val = Foam::zero();
    for(vector const& v : *this)
        val += v;
    val /= this->size();
    Info<<"| avg value: "<<val<<Foam::nl;
    Info<<"| size: "<<this->size()<<Foam::nl;
    Info<<"---------------------------------------------"<<Foam::nl;
}

void Foam::icoAdjointVelocityWallBC::set_dJdp_Wall
(
    std::function<Field<scalar>(const icoAdjointVelocityWallBC&)> expr
)
{
    Info<<"Set dJdp_Wall"<<Foam::nl;
    dJdp_Wall = expr;
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
        icoAdjointVelocityWallBC
    );
}


// ************************************************************************* //

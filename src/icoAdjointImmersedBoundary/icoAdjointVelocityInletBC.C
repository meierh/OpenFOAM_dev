#include "icoAdjointVelocityInletBC.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::icoAdjointVelocityInletBC::icoAdjointVelocityInletBC
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:fixedValueFvPatchField<vector>(p, iF, dict, false)
{}

Foam::icoAdjointVelocityInletBC::icoAdjointVelocityInletBC
(
    const icoAdjointVelocityInletBC& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fieldMapper& mapper
)
:fixedValueFvPatchField<vector>(ptf, p, iF, mapper)
{}

Foam::icoAdjointVelocityInletBC::icoAdjointVelocityInletBC
(
    const icoAdjointVelocityInletBC& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:fixedValueFvPatchField<vector>(pivpvf, iF)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void Foam::icoAdjointVelocityInletBC::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    if(!dJdp_Inlet)
        FatalErrorInFunction<<"dJdp_Inlet not set"<<exit(FatalError);    
    scalarField dJdp = dJdp_Inlet(*this);
    vectorField n = patch().nf();
    vectorField::operator=(n*-1*dJdp);
    
    fixedValueFvPatchField<vector>::updateCoeffs(); // sets updated_ to true
    
    
    //Info<<"---------------------------------------------"<<Foam::nl;
    //Info<<"| icoAdjointVelocityInletBC::updateCoeffs done"<<Foam::nl;
    /*
    vector val = Foam::zero();
    for(vector const& v : *this)
        val += v;
    val /= this->size();
    Info<<"inlet bc adj_U avg value: "<<val<<Foam::nl;
    */
    //Info<<"| size: "<<this->size()<<Foam::nl;
    //Info<<"---------------------------------------------"<<Foam::nl;
}

void Foam::icoAdjointVelocityInletBC::set_dJdp_Inlet
(
    std::function<Field<scalar>(const icoAdjointVelocityInletBC&)> expr
)
{
    Info<<"Set dJdp_Inlet"<<Foam::nl;
    dJdp_Inlet = expr;
}

/*
void Foam::icoAdjointVelocityInletBC::write(Ostream& os) const
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
        icoAdjointVelocityInletBC
    );
}


// ************************************************************************* //

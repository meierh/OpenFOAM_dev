#include "icoAdjointVelocityOutletBC.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::icoAdjointVelocityOutletBC::icoAdjointVelocityOutletBC
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:robinFvVectorPatchField(p, iF, dict),
nu("nu",dimensionSet(0,2,-1,0,0,0,0),0)
{}

Foam::icoAdjointVelocityOutletBC::icoAdjointVelocityOutletBC
(
    const icoAdjointVelocityOutletBC& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fieldMapper& mapper
)
:robinFvVectorPatchField(ptf, p, iF, mapper),
nu("nu",dimensionSet(0,2,-1,0,0,0,0),0)
{}

Foam::icoAdjointVelocityOutletBC::icoAdjointVelocityOutletBC
(
    const icoAdjointVelocityOutletBC& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:robinFvVectorPatchField(pivpvf, iF),
nu("nu",dimensionSet(0,2,-1,0,0,0,0),0)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void Foam::icoAdjointVelocityOutletBC::updateCoeffs()
{
    Info<<"icoAdjointVelocityOutletBC::updateCoeffs"<<Foam::nl;
    if (updated())
    {
        return;
    }
    
    write_a();
    write_b();
    write_c();    
    robinFvPatchField<vector>::updateCoeffs();    
    
    const fvPatch& thisPatch = icoAdjointVelocityOutletBC::patch();
    const vectorField n = thisPatch.nf().ref();
    
    const scalarField adjU_n = (*this) & n;
    const vectorField adjU_tt = (*this) - adjU_n*n;
    vectorField::operator=(adjU_tt);
}

void Foam::icoAdjointVelocityOutletBC::set_dJdu_Outlet
(
    std::function<Field<vector>(const icoAdjointVelocityOutletBC&)> expr
)
{
    dJdu_Outlet = expr;
}

void Foam::icoAdjointVelocityOutletBC::set_nu
(
    dimensionedScalar nu
)
{
    nu_set = true;
    this->nu=nu;
}

void Foam::icoAdjointVelocityOutletBC::write_a()
{
    const fvPatchField<vector>& u = patch().lookupPatchField<volVectorField,vector>("U");
    const tmp<vectorField> n = patch().nf();
    a = u&n; // u_n
}

void Foam::icoAdjointVelocityOutletBC::write_b()
{
    b = Field<scalar>(this->size(),nu.value());
}

void Foam::icoAdjointVelocityOutletBC::write_c()
{
    if(!dJdu_Outlet)
        FatalErrorInFunction<<"dJdu_Outlet not set"<<exit(FatalError);
    c = dJdu_Outlet(*this);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        icoAdjointVelocityOutletBC
    );
}


// ************************************************************************* //

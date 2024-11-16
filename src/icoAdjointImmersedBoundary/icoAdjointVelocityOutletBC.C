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
:robinFvVectorPatchField(p, iF, dict)
{}

Foam::icoAdjointVelocityOutletBC::icoAdjointVelocityOutletBC
(
    const icoAdjointVelocityOutletBC& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fieldMapper& mapper
)
:robinFvVectorPatchField(ptf, p, iF, mapper)
{}

Foam::icoAdjointVelocityOutletBC::icoAdjointVelocityOutletBC
(
    const icoAdjointVelocityOutletBC& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:robinFvVectorPatchField(pivpvf, iF)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void Foam::icoAdjointVelocityOutletBC::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    write_a();
    write_b();
    write_c();
    
    robinFvPatchField<vector>::updateCoeffs();
    
    const fvPatch& thisPatch = icoAdjointVelocityOutletBC::patch();
    const tmp<vectorField> n = thisPatch.nf();
    
    const tmp<scalarField> adjU_n = (*this) & n;
    const tmp<vectorField> adjU_tt = (*this) - adjU_n*n;
    
    vectorField::operator=(adjU_tt);
}

void Foam::icoAdjointVelocityOutletBC::set_dJdu_Outlet
(
    std::function<Field<vector>(const icoAdjointVelocityOutletBC&)> expr
)
{
    Info<<"Set dJdu_uOutlet"<<Foam::nl;
    dJdu_Outlet = expr;
}

void Foam::icoAdjointVelocityOutletBC::write_a()
{
    const fvPatchField<vector>& u = patch().lookupPatchField<volVectorField,vector>("U");
    const tmp<vectorField> n = patch().nf();
    //a = u&n; // u_n
}

void Foam::icoAdjointVelocityOutletBC::write_b()
{
    //b = 1;
}

void Foam::icoAdjointVelocityOutletBC::write_c()
{
    //c = vector(0,0,0);
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

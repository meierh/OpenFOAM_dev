#include "icoAdjointTemperatureOutletBC.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::icoAdjointTemperatureOutletBC::icoAdjointTemperatureOutletBC
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:robinFvScalarPatchField(p, iF, dict)
{
    Info<<"icoAdjointTemperatureOutletBC-------------------------------------------------------------"<<Foam::nl;
}

Foam::icoAdjointTemperatureOutletBC::icoAdjointTemperatureOutletBC
(
    const icoAdjointTemperatureOutletBC& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:robinFvScalarPatchField(ptf, p, iF, mapper)
{
    Info<<"icoAdjointTemperatureOutletBC-------------------------------------------------------------"<<Foam::nl;
}


Foam::icoAdjointTemperatureOutletBC::icoAdjointTemperatureOutletBC
(
    const icoAdjointTemperatureOutletBC& pivpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:robinFvScalarPatchField(pivpvf, iF)
{
    Info<<"icoAdjointTemperatureOutletBC-------------------------------------------------------------"<<Foam::nl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void Foam::icoAdjointTemperatureOutletBC::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    write_a();
    write_b();
    write_c();
    
    /*
    robinFvPatchField<vector>::updateCoeffs();
    
    const fvPatch& thisPatch = icoAdjointVelocityOutletBC::patch();
    const tmp<vectorField> n = thisPatch.nf();
    
    const tmp<scalarField> adjU_n = (*this) & n;
    const tmp<vectorField> adjU_tt = (*this) - adjU_n*n;
    
    vectorField::operator=(adjU_tt);
    */
}


void Foam::icoAdjointTemperatureOutletBC::set_dJdT_Outlet
(
    std::function<Field<scalar>(const icoAdjointTemperatureOutletBC&)> expr
)
{
    Info<<"Set dJdT_Outlet"<<Foam::nl;
    dJdT_Outlet = expr;
}

void Foam::icoAdjointTemperatureOutletBC::write_a()
{
    const fvPatchField<vector>& u = patch().lookupPatchField<volVectorField,vector>("U");
    const tmp<vectorField> n = patch().nf();
    //a = u&n; // u_n
}

void Foam::icoAdjointTemperatureOutletBC::write_b()
{
    //b = 1;
}

void Foam::icoAdjointTemperatureOutletBC::write_c()
{
    //c = 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        icoAdjointTemperatureOutletBC
    );
}


// ************************************************************************* //

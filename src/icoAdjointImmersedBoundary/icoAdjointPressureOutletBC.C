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
:fixedValueFvPatchField<scalar>(p, iF, dict, false),
nu("nu",dimensionSet(0,2,-1,0,0,0,0),0)
{}

Foam::icoAdjointPressureOutletBC::icoAdjointPressureOutletBC
(
    const icoAdjointPressureOutletBC& ptf,
    const fvPatch& p,
    const DimensionedField<scalar,volMesh>& iF,
    const fieldMapper& mapper
)
:fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
nu("nu",dimensionSet(0,2,-1,0,0,0,0),0)
{}

Foam::icoAdjointPressureOutletBC::icoAdjointPressureOutletBC
(
    const icoAdjointPressureOutletBC& pivpvf,
    const DimensionedField<scalar,volMesh>& iF
)
:fixedValueFvPatchField<scalar>(pivpvf, iF),
nu("nu",dimensionSet(0,2,-1,0,0,0,0),0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void Foam::icoAdjointPressureOutletBC::updateCoeffs()
{
    if(fixedValueFvPatchField<scalar>::updated())
    {
        return;
    }
    
    if(!nu_set)
        FatalErrorInFunction<<"nu not set"<<exit(FatalError);    
    if(!dJdu_Outlet)
        FatalErrorInFunction<<"dJdu_Outlet not set"<<exit(FatalError);
    
    const fvPatch& thisPatch = patch();
    
    //u_n * adj_u_n
    const fvsPatchField<scalar>& u_n = thisPatch.lookupPatchField<surfaceScalarField, scalar>("phi");
    const fvsPatchField<scalar>& adj_u_n = thisPatch.lookupPatchField<surfaceScalarField, scalar>("adj_phi");
    Field<scalar> u_n_adj_u_n = u_n*adj_u_n;
    
    //u & adj_u
    const fvPatchField<vector>& u = thisPatch.lookupPatchField<volVectorField, vector>("U");
    const fvPatchField<vector>& adj_u = thisPatch.lookupPatchField<volVectorField, vector>("adj_U");
    Field<scalar> u_dot_adj_u = u&adj_u;    

    tmp<Field<vector>> gradU = u.snGrad();
    vectorField n = thisPatch.nf();
    Field<scalar> gradUn = gradU&n;
    Field<scalar> nu_gradUn = nu.value()*gradUn;

    Field<scalar> dJdun = dJdu_Outlet(*this)&n;
    
    if(temperatureUsed)
    {
        const fvPatchField<scalar>& T = thisPatch.lookupPatchField<volScalarField, scalar>("T");
        const fvPatchField<scalar>& adj_T = thisPatch.lookupPatchField<volScalarField, scalar>("adj_T");
        Field<scalar> T_adj_T = T*adj_T;
        Field<scalar> adj_p_patch = u_n_adj_u_n + u_dot_adj_u + nu_gradUn + dJdun+T_adj_T;
        scalarField::operator==(adj_p_patch);
    }
    else
    {
        Field<scalar> adj_p_patch = u_n_adj_u_n + u_dot_adj_u + nu_gradUn + dJdun;
        scalarField::operator==(adj_p_patch);
    }

    fixedValueFvPatchField<scalar>::updateCoeffs();
}

void Foam::icoAdjointPressureOutletBC::set_dJdu_Outlet
(
    std::function<Field<vector>(const icoAdjointPressureOutletBC&)> expr
)
{
    dJdu_Outlet = expr;
}

void Foam::icoAdjointPressureOutletBC::set_nu
(
    dimensionedScalar nu
)
{
    nu_set = true;
    this->nu=nu;
}

void Foam::icoAdjointPressureOutletBC::set_temperatureUsed
(
    bool temperatureUsed
)
{
    temperatureUsed_set = true;
    this->temperatureUsed=temperatureUsed;
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

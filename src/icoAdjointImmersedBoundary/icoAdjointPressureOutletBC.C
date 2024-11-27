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
    
    //u & adj_u
    const fvPatchField<vector>& u = thisPatch.lookupPatchField<volVectorField, vector>("U");
    const fvPatchField<vector>& adj_u = thisPatch.lookupPatchField<volVectorField, vector>("adj_U");
    Field<scalar> u_dot_adj_u = u&adj_u;
    
    //u_n * adj_u_n
    scalarField u_n = u & thisPatch.nf();
    scalarField adj_u_n = adj_u & thisPatch.nf();
    Field<scalar> u_n_adj_u_n = u_n*adj_u_n;
    
    // nu dadj_udn
    tmp<Field<vector>> gradU = u.snGrad();
    vectorField n = thisPatch.nf();
    Field<scalar> gradUn = gradU&n;
    Field<scalar> nu_gradUn = nu.value()*gradUn;

    // dJdun
    Field<vector> dJdu = dJdu_Outlet(*this);
    Field<scalar> dJdun = dJdu & n;
    
    if(temperatureUsed)
    {
        const fvPatchField<scalar>& T = thisPatch.lookupPatchField<volScalarField, scalar>("T");
        const fvPatchField<scalar>& adj_T = thisPatch.lookupPatchField<volScalarField, scalar>("adj_T");
        Field<scalar> T_adj_T = T*adj_T;
        Field<scalar> adj_p_patch = u_n_adj_u_n + u_dot_adj_u /*+ nu_gradUn*/ + dJdun + T_adj_T;
        scalarField::operator=(adj_p_patch);
    }
    else
    {
        Field<scalar> adj_p_patch = u_n_adj_u_n + u_dot_adj_u /*+ nu_gradUn*/ + dJdun;
        scalarField::operator=(adj_p_patch);
    }

    fixedValueFvPatchField<scalar>::updateCoeffs();
    
    /*
    Info<<"||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||u_n_adj_u_n:"<<u_n_adj_u_n.size()<<Foam::nl;
    Info<<"u_n_adj_u_n:"<<u_n_adj_u_n<<Foam::nl;
    Info<<"||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||u_dot_adj_u:"<<u_dot_adj_u.size()<<Foam::nl;
    Info<<"u_dot_adj_u:"<<u_dot_adj_u<<Foam::nl;
    Info<<"||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||dJdun:"<<dJdun.size()<<Foam::nl;
    Info<<dJdun<<Foam::nl;
    */
    
    Info<<"---------------------------------------------temperatureUsed:"<<temperatureUsed<<Foam::nl;

    scalar avg_u_n = Foam::zero();
    for(scalar const& v : u_n)
        avg_u_n += v;
    avg_u_n /= this->size();
    Info<<" avg_u_n: "<<avg_u_n<<Foam::nl;
    
    scalar avg_adj_u_n = Foam::zero();
    for(scalar const& v : adj_u_n)
        avg_adj_u_n += v;
    avg_adj_u_n /= this->size();
    Info<<" avg_adj_u_n: "<<avg_adj_u_n<<Foam::nl;
    
    vector avg_u = Foam::zero();
    for(vector const& v : u)
        avg_u += v;
    avg_u /= this->size();
    Info<<" avg_u: "<<avg_u<<Foam::nl;
    
    vector avg_adj_u = Foam::zero();
    for(vector const& v : adj_u)
        avg_adj_u += v;
    avg_adj_u /= this->size();
    Info<<" avg_adj_u: "<<avg_adj_u<<Foam::nl;
    
    vector avg_dJdu = Foam::zero();
    for(vector const& v : dJdu)
        avg_dJdu += v;
    avg_dJdu /= this->size();
    Info<<" avg_dJdu: "<<avg_dJdu<<Foam::nl;
    
    Info<<"---------------------------------------------"<<Foam::nl;
    Info<<"| icoAdjointPressureOutletBC::updateCoeffs done"<<Foam::nl;
    scalar val = Foam::zero();
    for(scalar const& v : *this)
        val += v;
    val /= this->size();
    Info<<"| avg value: "<<val<<Foam::nl;
    Info<<"| size: "<<this->size()<<Foam::nl;
    Info<<"---------------------------------------------"<<Foam::nl;
    //Info<<(*this)<<Foam::nl;
    //FatalErrorInFunction<<"Temp stop"<<exit(FatalError);
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

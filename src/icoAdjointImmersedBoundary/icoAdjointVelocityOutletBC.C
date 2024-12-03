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
: fixedValueFvPatchField<vector>(p, iF, dict),
nu("nu",dimensionSet(0,2,-1,0,0,0,0),0)
{}

Foam::icoAdjointVelocityOutletBC::icoAdjointVelocityOutletBC
(
    const icoAdjointVelocityOutletBC& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fieldMapper& mapper
)
: fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
nu("nu",dimensionSet(0,2,-1,0,0,0,0),0)
{}

Foam::icoAdjointVelocityOutletBC::icoAdjointVelocityOutletBC
(
    const icoAdjointVelocityOutletBC& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
: fixedValueFvPatchField<vector>(pivpvf, iF),
nu("nu",dimensionSet(0,2,-1,0,0,0,0),0)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void Foam::icoAdjointVelocityOutletBC::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    /*
    write_a();
    write_b();
    write_c();    
    robinFvPatchField<vector>::updateCoeffs();    
    */
    
    //const fvPatchField<scalar>& adj_phi = patch().lookupPatchField<surfaceScalarField,scalar>("adj_phi"); // Here not set
    const fvPatchField<vector>& u = patch().lookupPatchField<volVectorField,vector>("U");
    scalarField u_n(mag(patch().nf() & u));
    vectorField adj_U;
    patchInternalField(adj_U);
    vectorField adj_U_n = patch().nf()*(patch().nf() & adj_U);
    
    vectorField dJdu = dJdu_Outlet(*this);
    vectorField dJdu_n = (dJdu & patch().nf())*patch().nf();
    vectorField dJdu_t = dJdu-dJdu_n;
    
    scalar epsilon = 1e-10;
    vectorField adj_U_t = dJdu_t/(u_n+epsilon);
    
    vectorField::operator=(adj_U_t + adj_U_n);
    //vectorField::operator=(adj_U_t + (adj_phi*patch().Sf() / sqr(patch().magSf())));
    
    fixedValueFvPatchField<vector>::updateCoeffs(); // sets updated_ to true
    
    /*
    Info<<"---------------------------------------------"<<Foam::nl;

    vector avgU = Foam::zero();
    for(vector const& v : u)
        avgU += v;
    avgU /= this->size();
    Info<<" avgU: "<<avgU<<Foam::nl;
    
    scalar avg_u_n = Foam::zero();
    for(scalar const& v : u_n)
        avg_u_n += v;
    avg_u_n /= this->size();
    Info<<" avg_u_n: "<<avg_u_n<<Foam::nl;
    
    vector avg_adj_U_n = Foam::zero();
    for(vector const& v : adj_U_n)
        avg_adj_U_n += v;
    avg_adj_U_n /= this->size();
    Info<<" avg_adj_U_n: "<<avg_adj_U_n<<Foam::nl;
    
    vector avg_dJdu = Foam::zero();
    for(vector const& v : dJdu)
        avg_dJdu += v;
    avg_dJdu /= this->size();
    Info<<" avg_dJdu: "<<avg_dJdu<<Foam::nl;
    */
    
    //Info<<"---------------------------------------------"<<Foam::nl;
    //Info<<"| icoAdjointVelocityOutletBC::updateCoeffs done"<<Foam::nl;
    /*
    vector val = Foam::zero();
    for(vector const& v : *this)
        val += v;
    val /= this->size();
    Info<<"outlet bc adj_U avg value: "<<val<<Foam::nl;
    */
    //Info<<"| size: "<<this->size()<<Foam::nl;
    //Info<<"---------------------------------------------"<<Foam::nl;
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
    nu.value()=0;
    this->nu=nu;
}

/*
void Foam::icoAdjointVelocityOutletBC::write_a()
{
    const fvPatchField<vector>& u = patch().lookupPatchField<volVectorField,vector>("U");
    const tmp<vectorField> n = patch().nf();
    a = mag(u&n); // u_n
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
*/

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

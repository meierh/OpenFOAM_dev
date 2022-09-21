#include "FSIFluidBoundaryConditions.H"

Foam::FSIFluidBoundaryConditions::FSIFluidBoundaryConditions
(
    const dimensionedScalar& alpha,
    const volScalarField& T,
    const volScalarField& p,
    const volVectorField& U,
    const fvMesh& mesh,
    const dimensionedScalar nu,
    const word& IBpatchName
):
alpha(alpha),
T(T),
p(p),
U(U),
mesh(mesh),
IBpatchName(IBpatchName),
nu(nu)
{
    auto cutCellBound = this->mesh.Sf().boundaryField();
    const fvBoundaryMesh& bound = mesh.boundary();
    IBPatchID = bound.findPatchID(IBpatchName);
    Info<<"BoundaryField size:"<<cutCellBound[IBPatchID].size()<<endl;
    Info<<"IB Patch id:"<<IBPatchID<<endl;
    
    computeIBHeatFlux();
    computeIBForce();
}

void Foam::FSIFluidBoundaryConditions::computeIBHeatFlux()
{
    tmp<GeometricField<double,fvsPatchField,surfaceMesh>> gradT = fvc::snGrad(T,IBpatchName);
    GeometricField<double,fvsPatchField,surfaceMesh>& gTF = gradT.ref();
    const GeometricField<double,fvsPatchField,surfaceMesh>::Boundary& gradTBoundary = gTF.boundaryField();
    const fvsPatchField<double> ibgradT = gradTBoundary[IBPatchID];
    
    Info<<"gradT size:"<<gradTBoundary.size()<<endl; 
    Info<<"ibgradT size:"<<ibgradT.size()<<endl;
}

void Foam::FSIFluidBoundaryConditions::computeIBForce()
{
    tmp<GeometricField<Tensor<double>,fvPatchField,volMesh>> gU = fvc::grad(U);
    GeometricField<Tensor<double>,fvPatchField,volMesh>& gradU = gU.ref();
    GeometricField<SymmTensor<double>,fvPatchField,volMesh> totalStress = symm(-p*tensor::one + nu*(gradU + gradU.T()));
    
    const vectorField& Sfp = mesh.Sf().boundaryField()[IBPatchID];
    const scalarField& magSfp = mesh.magSf().boundaryField()[IBPatchID];
    const symmTensorField& totalStressIB = totalStress.boundaryField()[IBPatchID];
    
    Field<Vector<double>> ibWallForces = (-Sfp/magSfp) & totalStressIB;
    Info<<"ibWallForces size:"<<ibWallForces.size()<<endl;

}

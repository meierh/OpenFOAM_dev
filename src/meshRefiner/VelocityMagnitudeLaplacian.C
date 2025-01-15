#include "VelocityMagnitudeLaplacian.H"

Foam::VelocityMagnitudeLaplacian::VelocityMagnitudeLaplacian
(
    fvMesh& mesh,
    LineStructure& structure,
    volScalarField& doRefine,
    dictionary& dynamicMeshDict,
    volVectorField& velocity
):
VelocityRefiner(mesh,structure,doRefine,dynamicMeshDict,velocity)
{}

/*
void Foam::VelocityRefiner::fieldRefinement()
{
    
    dimensionSet dimensions = fieldRefineDemands.dimensions();
    Foam::dimensioned<Foam::scalar> val("fieldRefinement",dimensions,-0.25);
    fieldRefineDemands = val;
    

    tmp<VolField<Tensor<scalar>>> velocityGradient = fvc::grad(velocity);
    for(label cellInd=0; cellInd<velocity.size(); cellInd++)
    {

    }
}
*/

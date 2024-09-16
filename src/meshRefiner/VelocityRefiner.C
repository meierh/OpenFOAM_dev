#include "VelocityRefiner.H"

Foam::VelocityRefiner::VelocityRefiner
(
    volScalarField& doRefine,
    LineStructure& structure,
    fvMesh& mesh,
    const dictionary& dynamicMeshDict,
    volVectorField& velocity
):
MeshRefiner(mesh,structure,doRefine,dynamicMeshDict),
velocity(velocity)
{}

void Foam::VelocityRefiner::fieldRefinement()
{
    /*
    dimensionSet dimensions = fieldRefineDemands.dimensions();
    Foam::dimensioned<Foam::scalar> val("fieldRefinement",dimensions,-0.25);
    fieldRefineDemands = val;
    */

    tmp<VolField<Tensor<scalar>>> velocityGradient = fvc::grad(velocity);
    for(label cellInd=0; cellInd<velocity.size(); cellInd++)
    {

    }
}

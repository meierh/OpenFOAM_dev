#include "RMSVelocityCurvature.H"

Foam::RMSVelocityCurvature::RMSVelocityCurvature
(
    volScalarField& doRefine,
    LineStructure& structure,
    fvMesh& mesh,
    const dictionary& dynamicMeshDict,
    volVectorField& velocity
):
VelocityRefiner(mesh,structure,doRefine,dynamicMeshDict,velocity)
{}

void Foam::RMSVelocityCurvature::fieldRefinement()
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

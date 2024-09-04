#include "VelocityRefiner.H"

Foam::VelocityRefiner::VelocityRefiner
(
    volScalarField& doRefine,
    const LineStructure& structure,
    fvMesh& mesh,
    Time& runTime,
    volVectorField& velocity,
    scalar refineGradient,
    scalar unrefineGradient
):
MeshRefiner(doRefine,structure,mesh,runTime),
velocity(velocity),
refineGradient(refineGradient),
unrefineGradient(unrefineGradient)
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

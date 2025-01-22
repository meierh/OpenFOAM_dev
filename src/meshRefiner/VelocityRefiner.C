#include "VelocityRefiner.H"

Foam::VelocityRefiner::VelocityRefiner
(
    fvMesh& mesh,
    LineStructure& structure,
    volScalarField& doRefine,
    dictionary& dynamicMeshDict,
    volVectorField& velocity
):
MeshRefiner(mesh,structure,doRefine,dynamicMeshDict),
velocity(velocity)
{}

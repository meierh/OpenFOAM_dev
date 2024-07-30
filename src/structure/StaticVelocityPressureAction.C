#include "StaticVelocityPressureAction.H"

Foam::StaticVelocityPressureAction::StaticVelocityPressureAction
(
    dynamicRefineFvMesh& mesh,
    LineStructure& structure,
    volVectorField& input_U,
    volVectorField& output_Uf
):
VelocityPressureForceInteraction(mesh,structure,input_U,output_Uf)
{
}

vector Foam::StaticVelocityPressureAction::getVelocity
(
    const LagrangianMarker* marker
)
{
    return vector(0,0,0);
}

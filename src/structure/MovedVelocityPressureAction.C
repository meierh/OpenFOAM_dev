#include "MovedVelocityPressureAction.H"

Foam::MovingVelocityPressureAction::MovingVelocityPressureAction
(
    const fvMesh& mesh,
    LineStructure& structure,
    const IOdictionary& stuctureDict,
    volVectorField& input_U,
    volVectorField& output_Uf
):
VelocityPressureForceInteraction(mesh,structure,input_U,output_Uf),
structureDict(structureDict),
movementDict(this->structureDict.subDict("rodMovement"))
{
}

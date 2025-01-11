#include "StaticVelocityPressureAction.H"

Foam::StaticVelocityPressureAction::StaticVelocityPressureAction
(
    const fvMesh& mesh,
    LineStructure& structure,
    volVectorField& input_U,
    volVectorField& output_Uf,
    const IOdictionary& structureDict,
    std::shared_ptr<MeshRefiner> refinement_,
    markerMeshType modusFieldToMarker,
    markerMeshType modusMarkerToField
):
VelocityPressureForceInteraction(mesh,structure,input_U,output_Uf,structureDict,refinement_,modusFieldToMarker,modusMarkerToField)
{}

void Foam::StaticVelocityPressureAction::preSolveMarkerMeshAdaption()
{
    if(!initalMarkerMeshAdaptionDone)
    {
        meshMarkerAdaptation();
        initalMarkerMeshAdaptionDone = true;
    }
}

bool Foam::StaticVelocityPressureAction::reconstructInterior
(
    bool firstIteration,
    scalar time,
    bool finalIteration
)
{
    return false;
}

Foam::vector Foam::StaticVelocityPressureAction::getCellVelocity
(
    label rodInd,
    scalar para,
    scalar angle,
    scalar radiusFrac
)
{
    return vector(0,0,0);
}

Foam::vector Foam::StaticVelocityPressureAction::getVelocity
(
    const LagrangianMarker* marker
)
{
    return vector(0,0,0);
}

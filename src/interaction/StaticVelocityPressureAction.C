#include "StaticVelocityPressureAction.H"

Foam::StaticVelocityPressureAction::StaticVelocityPressureAction
(
    const fvMesh& mesh,
    LineStructure& structure,
    volVectorField& input_U,
    volVectorField& output_Uf,
    std::shared_ptr<MeshRefiner> refinement_,
    markerMeshType modusFieldToMarker,
    markerMeshType modusMarkerToField
):
VelocityPressureForceInteraction(mesh,structure,input_U,output_Uf,refinement_,modusFieldToMarker,modusMarkerToField)
{}

void Foam::StaticVelocityPressureAction::preSolveMarkerMeshAdaption()
{
    if(!initalMarkerMeshAdaptionDone)
    {
        meshMarkerAdaptation();
        initalMarkerMeshAdaptionDone = true;
    }
    FatalErrorInFunction<<"Temp Stop"<<exit(FatalError);
}

Foam::vector Foam::StaticVelocityPressureAction::getVelocity
(
    const LagrangianMarker* marker
)
{
    return vector(0,0,0);
}

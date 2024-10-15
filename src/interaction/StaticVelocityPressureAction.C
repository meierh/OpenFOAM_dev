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
{
    Info<<" Create StaticVelocityPressureAction "<<Foam::endl;
    Info<<" Refinement set:"<<static_cast<bool>(refinement_)<<Foam::endl;
    if(refinement_)
    {
        Info<<"Do refinement"<<Foam::endl;
        refinement_->refineMeshOnStaticMarkers();
        refinement_->refineMeshAndMarkers();
    }
    structure.finalizeMarkers();
}

Foam::vector Foam::StaticVelocityPressureAction::getVelocity
(
    const LagrangianMarker* marker
)
{
    return vector(0,0,0);
}

#include "StaticVelocityPressureAction.H"

Foam::StaticVelocityPressureAction::StaticVelocityPressureAction
(
    const fvMesh& mesh,
    LineStructure& structure,
    volVectorField& input_U,
    volVectorField& output_Uf,
    markerMeshType modusFieldToMarker,
    markerMeshType modusMarkerToField
):
VelocityPressureForceInteraction(mesh,structure,input_U,output_Uf,nullptr,modusFieldToMarker,modusMarkerToField)
{
    Info<<" Create StaticVelocityPressureAction "<<Foam::endl;
    if(refinement_)
    {
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

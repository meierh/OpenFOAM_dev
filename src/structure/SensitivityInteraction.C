#include "SensitivityInteraction.H"

Foam::SensitivityInteraction::SensitivityInteraction
(
    dynamicRefineFvMesh& mesh,
    LineStructure& structure,
    markerMeshType modusFieldToMarker,
    markerMeshType modusMarkerToField
):
FieldMarkerStructureInteraction(mesh,structure,modusFieldToMarker,modusMarkerToField)
{
    Info<<"Completed SensitivityInteraction setup"<<Foam::endl;
}

#include "SensitivityInteraction.H"

Foam::SensitivityInteraction::SensitivityInteraction
(
    const fvMesh& mesh,
    LineStructure& structure,
    const IOdictionary& structureDict,
    markerMeshType modusFieldToMarker,
    markerMeshType modusMarkerToField
):
FieldMarkerStructureInteraction(mesh,structure,structureDict,modusFieldToMarker,modusMarkerToField)
{
    Info<<"Completed SensitivityInteraction setup"<<Foam::endl;
}

#include "SensitivityTemperatureInteraction.H"

Foam::SensitivityTemperatureInteraction::SensitivityTemperatureInteraction
(
    const fvMesh& mesh,
    LineStructure& structure,
    markerMeshType modusFieldToMarker,
    markerMeshType modusMarkerToField
):
SensitivityInteraction(mesh,structure,modusFieldToMarker,modusMarkerToField)
{}

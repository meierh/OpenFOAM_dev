#include "SensitivityVelocityPressureForceMarkerInteraction.H"

Foam::SensitivityVelocityPressureForceInteraction::SensitivityVelocityPressureForceInteraction
(
    const fvMesh& mesh,
    LineStructure& structure,
    markerMeshType modusFieldToMarker,
    markerMeshType modusMarkerToField
):
SensitivityInteraction(mesh,structure,modusFieldToMarker,modusMarkerToField)
{}

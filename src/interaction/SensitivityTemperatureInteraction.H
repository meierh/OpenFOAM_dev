#ifndef SENSITIVITYTEMPERATUREINTERACTION_H
#define SENSITIVITYTEMPERATUREINTERACTION_H
#include "SensitivityInteraction.H"

namespace Foam
{

class SensitivityTemperatureInteraction : public SensitivityInteraction
{
    public:
        explicit SensitivityTemperatureInteraction
        (
            const fvMesh& mesh,
            LineStructure& structure,
            markerMeshType modusFieldToMarker = markerMeshType::NonUniform,
            markerMeshType modusMarkerToField = markerMeshType::NonUniform
        );
        virtual ~SensitivityTemperatureInteraction()=default;
};
}
#endif

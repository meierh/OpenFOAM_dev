#ifndef SENSITIVITYVELOCITYPRESSUREINTERACTION_H
#define SENSITIVITYVELOCITYPRESSUREINTERACTION_H
#include "SensitivityInteraction.H"

namespace Foam
{

class SensitivityVelocityPressureForceInteraction : public SensitivityInteraction
{
    public:
        explicit SensitivityVelocityPressureForceInteraction
        (
            const fvMesh& mesh,
            LineStructure& structure,
            markerMeshType modusFieldToMarker = markerMeshType::NonUniform,
            markerMeshType modusMarkerToField = markerMeshType::NonUniform
        );
        virtual ~SensitivityVelocityPressureForceInteraction()=default;
};
}
#endif

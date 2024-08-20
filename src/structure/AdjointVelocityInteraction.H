#ifndef ADJOINTVELOCITYINTERACTION_H
#define ADJOINTVELOCITYINTERACTION_H
#include "SensitivityInteraction.H"

namespace Foam
{

class AdjointVelocityInteraction : public SensitivityInteraction
{
    public:
        explicit SensitivityInteraction
        (
            dynamicRefineFvMesh& mesh,
            LineStructure& structure,
            markerMeshType modusFieldToMarker = markerMeshType::NonUniform,
            markerMeshType modusMarkerToField = markerMeshType::NonUniform
        );

    protected:

};
#include "SensitivityInteraction.T.H"
}
#endif

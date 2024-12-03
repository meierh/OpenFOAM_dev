#ifndef SENSITIVITYVELOCITYPRESSURE_H
#define SENSITIVITYVELOCITYPRESSURE_H
#include "VelocityPressureForceMarkerInteraction.H"
#include "SensitivityInteraction.H"

namespace Foam
{

class SensitivityVelocityPressure : public SensitivityInteraction
{
    public:
        explicit SensitivityVelocityPressure
        (
            const fvMesh& mesh,
            LineStructure& structure,
            const VelocityPressureForceInteraction& primalInteraction,
            const IOdictionary& structureDict,
            markerMeshType modusFieldToMarker = markerMeshType::NonUniform,
            markerMeshType modusMarkerToField = markerMeshType::NonUniform
        );
        virtual ~SensitivityVelocityPressure()=default;
        
        scalar computeSensitivity(const Parameter& para) override;
        
    protected:
        const VelocityPressureForceInteraction& primalInteraction;
        volVectorField forcingDerivativeField;
        vector integrateVelocityForcingSensitivity(const Parameter& par);
        vector integrateVelocitySensitivity(const Parameter& par);
};
}
#endif

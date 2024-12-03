#ifndef SENSITIVITYTEMPERATURE_H
#define SENSITIVITYTEMPERATURE_H
#include "TemperatureInteraction.H"
#include "SensitivityInteraction.H"

namespace Foam
{

class SensitivityTemperature : public SensitivityInteraction
{
    public:
        explicit SensitivityTemperatureInteraction
        (
            const fvMesh& mesh,
            LineStructure& structure,
            const TemperatureInteraction& primalInteraction,
            const IOdictionary& structureDict,
            markerMeshType modusFieldToMarker = markerMeshType::NonUniform,
            markerMeshType modusMarkerToField = markerMeshType::NonUniform
        );
        virtual ~SensitivityTemperatureInteraction()=default;
        
        scalar computeSensitivity(const Parameter& para) override;
        
    protected:
        const TemperatureInteraction& primalInteraction;
        volScalarField heatingDerivativeField;
        scalar integrateTemperatureForcingSensitivity(const Parameter& par);
        scalar integrateTemperatureSensitivity(const Parameter& par);
};
}
#endif

#ifndef SENSITIVITYTEMPERATUREINTERACTION_H
#define SENSITIVITYTEMPERATUREINTERACTION_H
#include "TemperatureInteraction.H"
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
            const TemperatureInteraction& primalInteraction,
            volScalarField& input_adj_T,
            volScalarField& output_adj_Tf,
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
        
    private:
        volScalarField& input_adj_T;
        volScalarField& output_adj_Tf;
        
        DynamicList<scalar> markerFluidAdjointTemperature;
        DynamicList<scalar> makerCouplingAdjointHeating;
};
}
#endif

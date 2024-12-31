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
            volScalarField& adj_T,
            volScalarField& adj_fT,
            const IOdictionary& structureDict,
            markerMeshType modusFieldToMarker = markerMeshType::NonUniform,
            markerMeshType modusMarkerToField = markerMeshType::NonUniform
        );
        virtual ~SensitivityTemperatureInteraction()=default;
        
        scalar computeSensitivity(const Parameter& para) override;
        
        void solve(scalar timeStep);
        void solve(const pimpleSingleRegionControl&) override {};
        
        void store() override {FatalErrorInFunction<<"Not Implmented"<<exit(FatalError);}
        void setToTime(scalar time) override {FatalErrorInFunction<<"Not Implmented"<<exit(FatalError);}
        
        void recomputeMarkerValues() override;
        
    protected:
        const TemperatureInteraction& primalInteraction;
        volScalarField heatingDerivativeField;
        
        void interpolateAdjTemperatureToMarkers();
        void computeAdjCouplingHeatingOnMarkers(scalar timeStep);
        void interpolateAdjHeatingField();
        
        scalar integrateTemperatureForcingSensitivity(const Parameter& par);
        scalar integrateTemperatureSensitivity(const Parameter& par);
        
    private:
        volScalarField& adj_T;
        volScalarField& adj_fT;
        
        DynamicList<scalar> markerFluidAdjointTemperature;
        DynamicList<scalar> makerCouplingAdjointHeating;
};
}
#endif

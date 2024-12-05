#ifndef TEMPERATUREINTERACTION_H
#define TEMPERATUREINTERACTION_H
#include "FieldMarkerStructureInteraction.H"

namespace Foam
{
    
class TemperatureInteraction : public FieldMarkerStructureInteraction
{
    public:
        explicit TemperatureInteraction
        (
            const fvMesh& mesh,
            LineStructure& structure,
            volScalarField& input_T,
            volScalarField& output_Tf,
            const IOdictionary& structureDict,
            markerMeshType modusFieldToMarker = markerMeshType::NonUniform,
            markerMeshType modusMarkerToField = markerMeshType::NonUniform
        );
        
        void solve() override;
        
        void store() override;
        void setToTime(scalar time) override;
        
        const volScalarField& getTemperatureField() const {return input_T;}
        const volScalarField& getTemperatureForcingField() const {return output_Tf;}
        const DynamicList<scalar>& getMarkerCouplingHeating() const {return makerCouplingHeating;}
        
        void recomputeMarkerValues();

    protected:
        void interpolateTemperatureToMarkers();
        void computeCouplingHeatingOnMarkers();
        void computeRodHeating();
        void interpolateHeatingField();
        
        scalar sumHeating
        (
            std::function<bool(LagrangianMarker)> condition = 
                                [](LagrangianMarker marker){return true;}
        );
        
        scalar getTemperature(const LagrangianMarker* marker);

    private:
        volScalarField& input_T;
        volScalarField& output_Tf;
        
        DynamicList<scalar> markerFluidTemperature;
        DynamicList<scalar> makerCouplingHeating;
        DynamicList<scalar> rodHeating;
        
        std::map<scalar,std::tuple<DynamicList<scalar>,DynamicList<scalar>,DynamicList<scalar>>> storage;
        
        void assignForceOnRod();
};

}
#endif

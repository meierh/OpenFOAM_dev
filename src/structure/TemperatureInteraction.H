#ifndef TEMPERATUREINTERACTION_H
#define TEMPERATUREINTERACTION_H
#include "LineStructure.H"

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
            const IOdictionary& structureDict
        );
        
        void solve() override;
        
        void store() override;
        
        void setToTime(scalar time) override;

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
    
        const IOdictionary& structureDict;
        
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

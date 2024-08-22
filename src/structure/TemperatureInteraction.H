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
            fvMesh& mesh,
            LineStructure& structure,
            volScalarField& input_T,
            volScalarField& output_Tf
        );

        void interpolateTemperatureToMarkers();
        void computeCouplingHeatingOnMarkers();
        void computeRodHeating();
        void interpolateHeatingField();
        
        scalar sumHeating
        (
            std::function<bool(LagrangianMarker)> condition = 
                                [](LagrangianMarker marker){return true;}
        );
        
    protected:
        scalar getTemperature(const LagrangianMarker* marker);
        
    private:
        volScalarField& input_T;
        volScalarField& output_Tf;
        
        DynamicList<scalar> markerFluidTemperature;
        DynamicList<scalar> makerCouplingHeating;
        DynamicList<scalar> rodHeating;
        
        void assignForceOnRod();
};

}
#endif

#ifndef FIXEDTEMPERATUREACTION_H
#define FIXEDTEMPERATUREACTION_H
#include "TemperatureInteraction.H"

namespace Foam
{
    
class FixedTemperatureAction : public TemperatureInteraction
{
    public:
        explicit FixedTemperatureAction
        (
            const fvMesh& mesh,
            LineStructure& structure,
            volScalarField& input_T,
            volScalarField& output_Tf,
            const IOdictionary& structureDict,
            markerMeshType modusFieldToMarker = markerMeshType::NonUniform,
            markerMeshType modusMarkerToField = markerMeshType::NonUniform
        );
        
    protected:
        scalar getTemperature(const LagrangianMarker* marker);
        
    private:
        scalar readTemperatureFromDict();
        scalar fixedTemperature;
};

}
#endif

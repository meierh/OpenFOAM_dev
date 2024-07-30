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
            dynamicRefineFvMesh& mesh,
            LineStructure& structure,
            volScalarField& input_T,
            volScalarField& output_Tf,
            scalar temperature
        );
        
    protected:
        scalar getTemperature(const LagrangianMarker* marker);
        
    private:
        scalar fixedTemperature;
};

}
#endif

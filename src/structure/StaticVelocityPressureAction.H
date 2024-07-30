#ifndef STATICVELOCITYPRESSUREACTION_H
#define STATICVELOCITYPRESSUREACTION_H
#include "VelocityPressureForceMarkerInteraction.H"

namespace Foam
{

class StaticVelocityPressureAction : public VelocityPressureForceInteraction
{
    public:
        explicit StaticVelocityPressureAction
        (
            dynamicRefineFvMesh& mesh,
            LineStructure& structure,
            volVectorField& input_U,
            volVectorField& output_Uf
        );
        
    private:
        vector getVelocity(const LagrangianMarker* marker);

};

}
#endif

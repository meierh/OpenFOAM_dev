#ifndef FORCEDMOVEMENTVELOCITYPRESSUREACTION_H
#define FORCEDMOVEMENTVELOCITYPRESSUREACTION_H
#include "VelocityPressureForceMarkerInteraction.H"

namespace Foam
{

class ForcedMovementVelocityPressureAction : public VelocityPressureForceInteraction
{
    public:
        explicit ForcedMovementVelocityPressureAction
        (
            const fvMesh& mesh,
            LineStructure& structure,
            volVectorField& input_U,
            volVectorField& output_Uf
        );
        
    private:
        vector getVelocity(const LagrangianMarker* marker);

};

}
#endif

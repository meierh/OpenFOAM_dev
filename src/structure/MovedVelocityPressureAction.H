#ifndef MOVEDVELOCITYPRESSUREACTION_H
#define MOVEDVELOCITYPRESSUREACTION_H
#include "VelocityPressureForceMarkerInteraction.H"

namespace Foam
{

class MovingVelocityPressureAction : public VelocityPressureForceInteraction
{
    public:
        explicit MovedVelocityPressureAction
        (
            dynamicRefineFvMesh& mesh,
            LineStructure& structure,
            volVectorField& input_U,
            volVectorField& output_Uf
        );
};

#endif

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
            const fvMesh& mesh,
            LineStructure& structure,
            volVectorField& input_U,
            volVectorField& output_Uf
        );
        
        void preSolveMovement() override {}
        
    protected:
        std::unique_ptr<List<List<vector>>> getDeformation() override {return nullptr;}
        
    private:
        vector getVelocity(const LagrangianMarker* marker);

};

}
#endif

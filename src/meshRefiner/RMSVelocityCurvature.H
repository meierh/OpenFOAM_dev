#ifndef RMSVelocityCurvature_H
#define RMSVelocityCurvature_H
#include "VelocityRefiner.H"

namespace Foam
{

class RMSVelocityCurvature : public VelocityRefiner
{
    public:
        explicit RMSVelocityCurvature
        (
            fvMesh& mesh,
            LineStructure& structure,
            volScalarField& doRefine,
            dictionary& dynamicMeshDict,
            volVectorField& velocity
        );

    protected:
        //void fieldRefinement() override;
};

}
#endif

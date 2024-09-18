#ifndef RMSVELOCITYCURVATURE_H
#define RMSVELOCITYCURVATURE_H
#include "VelocityRefiner.H"

namespace Foam
{

class RMSVelocityCurvature : public VelocityRefiner
{
    public:
        explicit RMSVelocityCurvature
        (
            volScalarField& doRefine,
            LineStructure& structure,
            fvMesh& mesh,
            const dictionary& dynamicMeshDict,
            volVectorField& velocity
        );

    protected:
        void fieldRefinement() override;
};

}
#endif

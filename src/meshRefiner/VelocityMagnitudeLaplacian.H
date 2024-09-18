#ifndef VELOCITYREFINER_H
#define VELOCITYREFINER_H
#include "VelocityRefiner.H"
#include "fvcGrad.H"

namespace Foam
{

class VelocityMagnitudeLaplacian : public VelocityRefiner
{
    public:
        explicit VelocityMagnitudeLaplacian
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

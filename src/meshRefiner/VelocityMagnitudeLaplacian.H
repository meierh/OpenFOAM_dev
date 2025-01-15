#ifndef VelocityMagnitudeLaplacian_H
#define VelocityMagnitudeLaplacian_H
#include "VelocityRefiner.H"
#include "fvcGrad.H"

namespace Foam
{

class VelocityMagnitudeLaplacian : public VelocityRefiner
{
    public:
        explicit VelocityMagnitudeLaplacian
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

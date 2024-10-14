#ifndef VELOCITYREFINER_H
#define VELOCITYREFINER_H
#include "MeshRefiner.H"
#include "fvcGrad.H"

namespace Foam
{

class VelocityRefiner : public MeshRefiner
{
    public:
        explicit VelocityRefiner
        (
            volScalarField& doRefine,
            LineStructure& structure,
            fvMesh& mesh,
            dictionary& dynamicMeshDict,
            volVectorField& velocity
        );

    protected:
        scalar refineGradient;
        scalar unrefineGradient;
        volVectorField& velocity;
        void fieldRefinement() override;
};

}
#endif

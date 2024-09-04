#ifndef VELOCITYREFINER_H
#define VELOCITYREFINER_H
#include "MeshRefiner.H"
#include "fvcGrad.H"

namespace Foam
{

class VelocityRefiner : public MeshRefiner
{
    public:
        explicit VelocityRefiner(volScalarField& doRefine, const LineStructure& structure, fvMesh& mesh, Time& runTime, volVectorField& velocity);

    protected:
        scalar refineGradient;
        scalar unrefineGradient;
        volVectorField& velocity;
        void fieldRefinement() override;
};

}
#endif

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
            fvMesh& mesh,
            LineStructure& structure,
            volScalarField& doRefine,
            dictionary& dynamicMeshDict,
            volVectorField& velocity
        );
        
    protected:
        scalar refineGradient;
        scalar unrefineGradient;
        volVectorField& velocity;
};

}
#endif

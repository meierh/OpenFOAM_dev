#ifndef VelocityJump_H
#define VelocityJump_H
#include "VelocityRefiner.H"

namespace Foam
{
class VelocityJump : public VelocityRefiner
{
    public:
        explicit VelocityJump
        (
            fvMesh& mesh,
            LineStructure& structure,
            volScalarField& doRefine,
            dictionary& dynamicMeshDict,
            volVectorField& velocity
        );
        
    protected:
        void fieldRefinement() override;
        scalar refineJumpLimit;
        scalar unrefineJumpLimit;
        scalar timeOffset;
        bool timeOffsetSet;
        label maxRefinement;
        bool maxRefinementSet;
        volScalarField maxCellJump;
};
}
#endif

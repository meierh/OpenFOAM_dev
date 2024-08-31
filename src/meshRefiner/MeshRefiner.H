#ifndef MESHREFINER_H
#define MESHREFINER_H
#include "fvMesh.H"
#include "LineStructure.H"

namespace Foam
{

class MeshRefiner
{
    public:
        explicit MeshRefiner(volScalarField& doRefine, const LineStructure& structure, fvMesh& mesh, Time& runTime);
        virtual ~MeshRefiner() = default;

        virtual void adaptMesh();
    protected:
        const scalar REFINE   = +1;
        const scalar MUSTKEEP = 0.25;
        const scalar DONTCARE = -0.25;
        const scalar UNREFINE = -1;

        const scalar REFINE_LIM   = +0.5;
        const scalar MUSTKEEP_LIM =  0;
        const scalar DONTCARE_LIM = -0.5;

        fvMesh& mesh;
        Time& runTime;
        const LineStructure& structure;
        volScalarField& doRefine;

        volScalarField fieldRefineDemands;
        virtual void fieldRefinement();

        volScalarField markerRefineDemands;
        void markerRefinement(scalar defaultValue);

    private:
        scalar refinementDemandMerge(scalar fieldDemand, scalar markerDemand);
};

}
#endif

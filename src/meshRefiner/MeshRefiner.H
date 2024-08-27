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
        volScalarField& doRefine;
        const LineStructure& structure;
        fvMesh& mesh;
        Time& runTime;

        void setMarkerRefinementInfo();
        volScalarField markerRefinementInfo;
};

}
#endif

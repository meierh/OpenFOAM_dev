#ifndef MESHREFINER_H
#define MESHREFINER_H
#include "fvCFD.H"
#include "LineStructure.H"

namespace Foam
{
    
class MeshRefiner
{
    public:
        explicit MeshRefiner(volScalarField& doRefine, const LineStructure& structure, fvMesh& mesh);
        virtual ~MeshRefiner() = default;
        
    protected:
        volScalarField& doRefine;
        const LineStructure& structure;
        fvMesh& mesh;
        
        void setRefinementCells();
};

}
#endif

#ifndef FVMESHOCTREE_H
#define FVMESHOCTREE_H
#include <memory>
#include "fvCFD.H"

namespace Foam
{
class OctreeNode
{
    public:
        vector seperationCoord;
        std::unique_ptr<OctreeNode> childs[2][2][2];
        DynamicList<label> cells;
};
    
class FvMeshOctree
{
    public:
        explicit FvMeshOctree(const fvMesh& mesh);
        FvMeshOctree(const FvMeshOctree&) = delete;
                
    private:
        const fvMesh& mesh;
        void computingBoundingBox();
        void constructOctree(const fvMesh& mesh);
        void createNode(const fvMesh& mesh,std::unique_ptr<OctreeNode>& node,uint maxNumberOfCellsInNode);
        
        vector boundingBoxTop;
        vector boundingBoxBottom;
        
        std::unique_ptr<OctreeNode> root;
        std::pair<vector,vector> topBottomFromPoints(const pointField& points) const;
        std::pair<vector,vector> topBottomFromPoints(const DynamicList<vector>& points) const;
        void assignCellToSide
        (
            label cellInd,
            vector splitVec,
            DynamicList<std::tuple<label,label,label>>& bins
        );
        std::tuple<DynamicList<label>,DynamicList<label>,DynamicList<label>> assignPointsToSide(const pointField& points,vector splitVec)const;
        
};
}
#endif

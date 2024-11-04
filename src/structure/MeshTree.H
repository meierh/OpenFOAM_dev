#ifndef MESHTREE_H
#define MESHTREE_H

#include <memory>
#include <unordered_set>
#include "fvMesh.H"
#include "activeRodMesh.h"
#include "BoundingBox.H"

namespace Foam
{
    
class MeshTree
{
    public:
        MeshTree(const fvMesh& mesh);
        
    private:
        const fvMesh& mesh;

        struct node;
        using CellData = std::tuple<label,BoundingBox,std::vector<node*>>;
        template <typename T>
        using OctoSplit = FixedList<FixedList<FixedList<T,2>,2>,2>;
        struct node
        {
            vector split;
            BoundingBox box;
            OctoSplit<std::unique_ptr<node>> leafs;
            DynamicList<CellData*> data;
        };
        
        std::unique_ptr<node> root;
        DynamicList<CellData> cellBBs;
        
        void computeCellBBs();
        void constructTree();
        BoundingBox minMaxSpan(const DynamicList<CellData*>& boxes) const;
        FixedList<FixedList<FixedList<bool,2>,2>,2> splitOutNode(node* oneNode);
        labelList side(const CellData* cData, scalar value, label dim) const;
};

}

#endif

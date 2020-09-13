#ifndef KDTREE_H
#define KDTREE_H

#include "fvCFD.H"
#include <memory>
#include <set>
#include "Nurbs.H"

namespace Foam
{   
class KdTree
{
    private:
        struct Node
        {
            Node* parent;
            Node* left;
            Node* right;
            BoundingBox MinMaxBox;
            labelList nurbsCurves;
            label divideDim;
            scalar divideBound;
        };
                
        Node* _nil;
        Node* root;
                
        const std::unique_ptr<List<Nurbs*>> Items;
        const label maxCurvesPerNode;

        List<Foam::BoundingBox> listMinMaxBoxes;
        
    public:
        explicit KdTree
        (
            std::unique_ptr<List<Nurbs*>> Items,
            label maxCurvesPerNode=1
        );
        ~KdTree();
        labelList nearNurbsCurves(vector point);
        
    private:
        Node *newNode(Node* parent);
        void constructTree
        (
            Node* thisNode,
            labelList nurbsCurves,
            label treeHeight,
            labelList firstLevel,
            labelList secondLevel,
            labelList thirdLevel
        );
        bool isInsideMinMaxBox
        (
            Foam::BoundingBox MinMaxBox,
            vector point
        );
        void traverseKdTree
        (
            Node* currentNode,
            vector point,
            std::set<label>* nearNurbsList
        ) const;
        void recursiveNodeDeleter(Node* thisNode);
        
        friend void UnitTest_KdTree();
};

void UnitTest_KdTree();
}
#endif

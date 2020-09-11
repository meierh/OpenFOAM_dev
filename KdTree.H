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
        
        label maxCurvesPerNode;
        
        Node* _nil;
        Node* root;
        
    private:
        std::unique_ptr<List<Nurbs*>> Items;
        List<Foam::BoundingBox> listMinMaxBoxes;
        
    public:
        explicit KdTree(std::unique_ptr<List<Nurbs*>> Items,label maxCurvesPerNode);
        labelList nearNurbsCurves(vector point);
        
    private:
        Node *newNode(Node* parent);
        void constructTree
        (
            Node* thisNode,
            labelList nurbsCurves,
            label treeHeight,
            labelList firstLevel,
            labelList secondLevel
        );
        bool isInsideMinMaxBox(Foam::BoundingBox MinMaxBox, vector point);
        void traverseKdTree(Node* currentNode,vector point,std::set<label>* nearNurbsList);
};

void UnitTest_KdTree();
}
#endif

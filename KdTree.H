#ifndef KDTREE_H
#define KDTREE_H

#include "fvCFD.H"
#include <memory>
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
                
        List<std::shared_ptr<Nurbs>> Items;
        label maxCurvesPerNode;

        List<Foam::BoundingBox> listMinMaxBoxes;
        
    public:
        explicit KdTree
        (
            List<std::shared_ptr<Nurbs>> Items,
            label maxCurvesPerNode=1
        );
        
        KdTree() {}
        
        ~KdTree();
        std::unique_ptr<labelList> nearNurbsCurves(vector point) const;
        
        labelList printPath(vector point);
        
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
        void traverseKdTree
        (
            Node* currentNode,
            vector& point,
            std::unique_ptr<labelList>& nearNurbsList
        ) const;
        
        void printTraversedKdTree
        (
            Node* currentNode,
            vector point,
            labelList& nearNurbsList
        ) const;
        
        void recursiveNodeDeleter(Node* thisNode);
        
        friend void UnitTest_KdTree();
        friend bool testTreeForDuplicatesOnPath
        (
            const KdTree::Node* thisNode,
            labelList foundCurves,
            const KdTree::Node* _nil
        );
        friend bool nonIncreasingBoundingBox(const KdTree::Node* thisNode, const KdTree::Node* _nil);
};

void UnitTest_KdTree();
bool testTreeForDuplicatesOnPath(const KdTree::Node* thisNode,labelList foundCurves, const KdTree::Node* _nil);
bool nonIncreasingBoundingBox(const KdTree::Node* thisNode, const KdTree::Node* _nil);
}
#endif

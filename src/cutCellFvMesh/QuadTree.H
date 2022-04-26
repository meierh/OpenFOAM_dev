#ifndef QUADTREE_H
#define QUADTREE_H

#include "fvCFD.H"
#include <memory>
#include "BsTree.H"
#include "Nurbs2D.H"

namespace Foam
{
class QuadTree : public BsTree
{
    private:
        struct Node
        {
            Node* parent;
            Node* leftU_leftV;
            Node* rightU_leftV;
            Node* leftU_rightV;
            Node* rightU_rightV;
            scalar minU;
            scalar maxU;
            scalar minV;
            scalar maxV;
            scalar divideBoundU;
            scalar divideBoundV;
            BoundingBox MinMaxBox;
        };
        const Nurbs2D Surface;
        
        Node* _nil;
        Node* root;
        
    public:
        QuadTree
        (
            Nurbs2D Surface,
            label nbrSplitsBetweenCPs = 1
        );
        //~QuadTree();
        
        DynamicList<FixedList<scalar,2>> nearestPoints(vector point) const;

        //scalar closestParaOnNurbsToPoint(vector point) const = delete;
        FixedList<scalar,2> closestParaOnNurbsToPoint(vector point) const;
                
    private:
        Node *newNode(Node* parent = nullptr,scalar min=0, scalar max=0)=delete;
        Node *newNode(Node* parent = nullptr,scalar minU=0, scalar maxU=0, scalar minV=0, scalar maxV=0);
        void constructTree(Node* thisNode, int height=0);
        void recursiveNodeDeleter(Node* thisNode);
               
        void traverseBsTree
        (
            Node* currentNode,
            vector point,
            DynamicList<FixedList<scalar,2>>& coordNurbs
        ) const;
        
        //friend void UnitTest_BsTree();
        friend class TESTQUADTREE;
};

//void UnitTest_BsTree();
}

#endif

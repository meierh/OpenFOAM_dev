#ifndef UNITTESTQUADTREE_H
#define UNITTESTQUADTREE_H

#include "QuadTree.H"
#include "UnitTestBsTree.H"
#include <math.h>
#include <cassert>
#include <random>
#include <vector>

namespace Foam
{
class TESTQUADTREE : public TESTBSTREE
{
    public:
        TESTQUADTREE()
        {
            nearestLeaf();
            Info<<"QuadTree Test done"<<endl;
        }
        
    protected:
        void nearestLeaf();
        void collectLeafCenter(DynamicList<FixedList<scalar,2>>& leafCenters, const QuadTree::Node* const node, const QuadTree& tree);
        void collectLeafCenterWithPointInBox(DynamicList<FixedList<scalar,2>>& leafCenters, const QuadTree::Node* const node, const QuadTree& tree, const vector point);
        bool checkNearestPoints(const QuadTree& tree, DynamicList<FixedList<scalar,2>> nearestPoints, vector point);
        bool checkNearestNeighbor(const QuadTree& tree, FixedList<scalar,2> checkNearestNeighbor, vector point);
};
}
#endif

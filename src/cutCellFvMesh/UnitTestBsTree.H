#ifndef UNITTESTBSTREE_H
#define UNITTESTBSTREE_H

#include "BsTree.H"
#include <math.h>
#include <cassert>
#include <random>

namespace Foam
{
class TESTBSTREE
{
    public:
        TESTBSTREE()
        {
            nearestLeaf();
            Info<<"BsTree Test done"<<endl;
        }
        
    protected:
        void nearestLeaf();
        void collectLeafCenter(DynamicList<scalar>& leafCenters, const BsTree::Node* const node, const BsTree& tree);
        bool checkNearestPoints(const BsTree& tree, scalarList nearestPoints, vector point);
        bool checkNearestNeighbor(const BsTree& tree, scalar checkNearestNeighbor, vector point);
};
}
#endif

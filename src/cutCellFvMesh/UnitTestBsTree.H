#ifndef UNITTESTNURBS2D_H
#define UNITTESTNURBS2D_H

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
};
}
#endif

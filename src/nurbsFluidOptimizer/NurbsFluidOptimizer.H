#ifndef NURBSFLUIDOPTIMIZER_H
#define NURBSFLUIDOPTIMIZER_H

#include <stdio.h>
#include <iostream>
//#include <direct.h>
#include <random>
#include <chrono>
#include <memory>
#include "dynamicRefineFvMesh.H"

//#include "activeRodMesh.h"
//#include "rodOpt.h"
//#include "rodTools.h"
//#include "rodLatticePBC.h"
//#include "rodImport.h"
//#include "rodCScircle.h"
//#include "rodCSrectangle.h"
//#include "timing.h"
//#include "rodMesh.h"
//#include "bvh.h"

//#include <dirent.h>

#define _USE_MATH_DEFINES
//#include <math.h>

#include "cutCellFvMesh.H"

namespace Foam
{
class NurbsFluidOptimizer
{
    public:
        explicit NurbsFluidOptimizer
        ( 
            Foam::argList args
        );
        ~NurbsFluidOptimizer();
        NurbsFluidOptimizer(const NurbsFluidOptimizer&) = delete;
        
        void updateNurbsParameters();
        std::shared_ptr<std::vector<NurbsData>> getNurbsParameters();
        bool loop();
        
    private:
        std::shared_ptr<std::vector<NurbsData>> nurbsGeom;
};
}
#endif

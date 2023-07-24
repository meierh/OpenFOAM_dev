#ifndef TFSITAPESTORAGE_H
#define TFSITAPESTORAGE_H

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
class TfsiTapeStorage
{
    public:
        explicit TfsiTapeStorage
        ( 
            cutCellFvMesh& mesh 
        );
        
        void registerField(volScalarField& field);
        void registerField(volVectorField& field);
        void record(scalar timeValue);
        
    private:
        cutCellFvMesh& mesh;
        std::stack<pointField> verticeTape;
        std::stack<scalar> times;
        
        template<typename T>
        class Tape
        {
            T& field;
            std::stack<T> fieldTape;
            std::stack<scalar> times;
            Tape(T& field):field(field){};
            void push(scalar time){fieldTape.push(field);times.push(time);};
            void pop(){fieldTape.push(field);times.push(time);};
        }
        std::vector<Tape<volScalarField>> scalarTape;
        std::vector<Tape<volVectorField>> vectorTape;
};
}
#endif

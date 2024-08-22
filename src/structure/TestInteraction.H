#ifndef TESTINTERACTION_H
#define TESTINTERACTION_H
#include "LineStructure.H"

namespace Foam
{

class TestInteraction : public FieldMarkerStructureInteraction
{
    public:
        explicit TestInteraction
        (
            fvMesh& mesh,
            volScalarField& testField
        );
        
        void printSupportToField();
        
    private:
        volScalarField& testField;

};

}
#endif

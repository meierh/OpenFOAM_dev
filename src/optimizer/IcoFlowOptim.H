#ifndef ICOFLOWOPTIM_H
#define ICOFLOWOPTIM_H
#include "IcoFlow.H"

namespace Foam
{
    
class IcoFlowOptim : public IcoFlow
{
    public:
        explicit IcoFlowOptim(int argc,char *argv[]);
        virtual ~IcoFlowOptim()=default;
        
    protected:
        // --- create Piso --- //
        void createPiso();

        Foam::pisoControl pisoAdjoint;
        // --- ----------- --- //
        
        // --- create Fields --- //
        void createFields() override;
        
        volScalarField adj_p;
        volVectorField adj_U;
        surfaceScalarField phi;
        label pRefCell;
        scalar pRefValue;
        // --- ----------- --- //
};

}
#endif

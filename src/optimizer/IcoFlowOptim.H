#ifndef ICOFLOWOPTIM_H
#define ICOFLOWOPTIM_H
#include "IcoFlow.H"
#include "SensitivityInteraction.H"

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
        surfaceScalarField adj_phi;
        label adj_pRefCell;
        scalar ajd_pRefValue;
        // --- ----------- --- //
};

}
#endif

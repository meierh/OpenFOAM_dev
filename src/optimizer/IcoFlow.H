#ifndef ICOFLOW_H
#define ICOFLOW_H
#include "Application.H"
#include "pisoControl.H"

namespace Foam
{
    
class IcoFlow : public Application
{
    public:
        explicit IcoFlow(int argc,char *argv[]);
        virtual ~IcoFlow()=default;
        
    protected:
        // --- create Piso --- //
        void createPiso();

        Foam::pisoControl piso;
        // --- ----------- --- //
        
        // --- create Fields --- //
        void createFields() override;
        
        IOdictionary transportProperties;
        dimensionedScalar nu;
        volScalarField p;
        volVectorField U;
        surfaceScalarField phi;
        label pRefCell;
        scalar pRefValue;
        // --- ----------- --- //
        
        // --- create Fields --- //
        void initContinuityErrs();
        
        scalar cumulativeContErr;
        // --- ----------- --- //
};

}
#endif

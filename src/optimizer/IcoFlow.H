#ifndef ICOFLOW_H
#define ICOFLOW_H
#include "Application.H"
#include "pisoControl.H"
#include <memory>

#include "VelocityPressureForceInteraction.H"
#include "TemperatureInteraction.H"

namespace Foam
{
    
class IcoFlow : public Application
{
    public:
        explicit IcoFlow(int argc,char *argv[]);
        virtual ~IcoFlow()=default;
        
        virtual void run();
        
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
        
        // --- create Immersed Boundary Fields --- //
        virtual void setupImmersedBoundary();
        std::unique_ptr<VelocityPressureForceInteraction> velocityImmersedBoundary;
        std::unique_ptr<volVectorField> f_U;
        std::unique_ptr<TemperatureInteraction> temperatureImmersedBoundary;
        std::unique_ptr<volVectorField> f_U;
        // --- ----------- --- //
};

}
#endif

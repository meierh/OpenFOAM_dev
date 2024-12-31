#ifndef SENSITIVITYSTATICVELOCITYPRESSUREACTION_H
#define SENSITIVITYSTATICVELOCITYPRESSUREACTION_H
#include "SensitivityVelocityPressureForceInteraction.H"

namespace Foam
{

class SensitivityVelocityPressureAction : public SensitivityVelocityPressureForceInteraction
{
    public:
        explicit SensitivityVelocityPressureAction
        (
            const fvMesh& mesh,
            LineStructure& structure,
            const VelocityPressureForceInteraction& primalInteraction,
            volVectorField& input_adj_U,
            volVectorField& output_adj_Uf,
            const IOdictionary& structureDict,
            markerMeshType modusFieldToMarker = markerMeshType::NonUniform,
            markerMeshType modusMarkerToField = markerMeshType::NonUniform
        );
        virtual ~SensitivityVelocityPressureAction()=default;
        
        scalar computeSensitivity(const Parameter& para) override;
        
    protected:
        const VelocityPressureForceInteraction& primalInteraction;
        volVectorField forcingDerivativeField;
        vector integrateVelocityForcingSensitivity(const Parameter& par);
        vector integrateVelocitySensitivity(const Parameter& par);
        
    private:
        volVectorField& input_adj_U;
        volVectorField& output_adj_Uf;
        
        DynamicList<vector> markerFluidAdjointVelocity;
        DynamicList<vector> makerCouplingAdjointForce;

};
}
#endif

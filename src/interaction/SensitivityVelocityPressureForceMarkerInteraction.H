#ifndef SENSITIVITYVELOCITYPRESSUREINTERACTION_H
#define SENSITIVITYVELOCITYPRESSUREINTERACTION_H
#include "VelocityPressureForceMarkerInteraction.H"
#include "SensitivityInteraction.H"

namespace Foam
{

class SensitivityVelocityPressureForceInteraction : public SensitivityInteraction
{
    public:
        explicit SensitivityVelocityPressureForceInteraction
        (
            const fvMesh& mesh,
            LineStructure& structure,
            const VelocityPressureForceInteraction& primalInteraction,
            volVectorField& adj_U,
            volVectorField& adj_fU,
            const IOdictionary& structureDict,
            markerMeshType modusFieldToMarker = markerMeshType::NonUniform,
            markerMeshType modusMarkerToField = markerMeshType::NonUniform
        );
        virtual ~SensitivityVelocityPressureForceInteraction()=default;
        
        scalar computeSensitivity(const Parameter& para) override;
        
        void solve(scalar timeStep);
        void solve(bool firstIteration,scalar time,bool finalIteration) override {};
        
        void store() override {FatalErrorInFunction<<"Not Implmented"<<exit(FatalError);}
        void setToTime(scalar time) override {FatalErrorInFunction<<"Not Implmented"<<exit(FatalError);}
        
        void recomputeMarkerValues() override;
        
    protected:
        const VelocityPressureForceInteraction& primalInteraction;
        volVectorField forcingDerivativeField;
        
        void interpolateAdjVelocityToMarkers();
        void computeAdjCouplingForceOnMarkers(scalar timeStep);
        void interpolateAdjFluidForceField();
        
        vector integrateVelocityForcingSensitivity(const Parameter& par);
        vector integrateVelocitySensitivity(const Parameter& par);
        
    private:
        volVectorField& adj_U;
        volVectorField& adj_fU;
        
        DynamicList<vector> markerFluidAdjointVelocity;
        DynamicList<vector> makerCouplingAdjointForce;        
};
}
#endif

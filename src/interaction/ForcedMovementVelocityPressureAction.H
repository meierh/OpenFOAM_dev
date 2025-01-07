#ifndef FORCEDMOVEMENTVELOCITYPRESSUREACTION_H
#define FORCEDMOVEMENTVELOCITYPRESSUREACTION_H
#include "VelocityPressureForceMarkerInteraction.H"
#include <dlfcn.h>

namespace Foam
{

class ForcedMovementVelocityPressureAction : public VelocityPressureForceInteraction
{
    public:
        explicit ForcedMovementVelocityPressureAction
        (
            const fvMesh& mesh,
            LineStructure& structure,
            volVectorField& input_U,
            volVectorField& output_Uf,
            const IOdictionary& structureDict,
            std::shared_ptr<MeshRefiner> refinement_,
            markerMeshType modusFieldToMarker = markerMeshType::NonUniform,
            markerMeshType modusMarkerToField = markerMeshType::NonUniform
        );
               
        void preSolveMovement() override;
        void preSolveMarkerMeshAdaption() override;
        
        void subTimestepStructureMovement() override {};
        void subTimestepMarkerMeshAdaption() override {};
        
        void recomputeMarkerValues() override {NotImplemented;}
        
    protected:
        std::unique_ptr<List<List<vector>>> readDeformationDict();
        
    private:
        vector getVelocity(const LagrangianMarker* marker) override;
        void constructMovementFunction();
        //std::unique_ptr<List<List<vector>>> readMovementFromDict();
};

}
#endif

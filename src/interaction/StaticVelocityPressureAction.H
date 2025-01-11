#ifndef STATICVELOCITYPRESSUREACTION_H
#define STATICVELOCITYPRESSUREACTION_H
#include "VelocityPressureForceMarkerInteraction.H"

namespace Foam
{

class StaticVelocityPressureAction : public VelocityPressureForceInteraction
{
    public:
        explicit StaticVelocityPressureAction
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
        
        void preSolveMovement() override {}
        void preSolveMarkerMeshAdaption() override;
        
        void subTimestepStructureMovement() override {};
        void subTimestepMarkerMeshAdaption() override {};
                
    protected:
        bool reconstructInterior(bool firstIteration,scalar time,bool finalIteration) override;
        vector getCellVelocity(label rodInd,scalar para,scalar angle,scalar radiusFrac) override;
        
    private:
        bool initalMarkerMeshAdaptionDone = false;
        vector getVelocity(const LagrangianMarker* marker);

};

}
#endif

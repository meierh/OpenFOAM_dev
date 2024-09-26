#ifndef VELOCITYPRESSUREFORCEMARKERINTERACTION_H
#define VELOCITYPRESSUREFORCEMARKERINTERACTION_H
#include "FieldMarkerStructureInteraction.H"
#include "MeshRefiner.H"

namespace Foam
{

class VelocityPressureForceInteraction : public FieldMarkerStructureInteraction
{
    public:
        explicit VelocityPressureForceInteraction
        (
            const fvMesh& mesh,
            LineStructure& structure,
            volVectorField& input_U,
            volVectorField& output_Uf,
            std::shared_ptr<MeshRefiner> refinement_,
            markerMeshType modusFieldToMarker = markerMeshType::NonUniform,
            markerMeshType modusMarkerToField = markerMeshType::NonUniform
        );

        virtual void preSolveMovement()=0;
        
        void solve() override;
        
        void store() override;
        void setToTime(scalar time) override;
        
        const volVectorField& getReferenceInOutField() const {return input_U;}
        const DynamicList<vector>& getMarkerCouplingForce() const {return markerCouplingForce;}
        
    protected:
        void interpolateFluidVelocityToMarkers();
        void computeCouplingForceOnMarkers();
        virtual void computeRodForceMoment();
        void moveMarkersAndAdaptMesh();
        void interpolateFluidForceField();
        
        vector sumForces
        (
            std::function<bool(LagrangianMarker)> condition = 
                                [](LagrangianMarker marker){return true;}
        );
        
        volVectorField& input_U;
        volVectorField& output_Uf;
        std::shared_ptr<MeshRefiner> refinement_;
        
        virtual vector getVelocity(const LagrangianMarker* marker) = 0;
        
        DynamicList<vector> markerFluidVelocity;
        DynamicList<vector> markerCouplingForce;
        DynamicList<vector> rodForce;
        DynamicList<vector> rodMoment;
        
        std::map<scalar,std::tuple<DynamicList<vector>,DynamicList<vector>,DynamicList<vector>,DynamicList<vector>>> storage;
};
}
#endif

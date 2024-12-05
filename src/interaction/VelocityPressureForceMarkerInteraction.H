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
            const IOdictionary& structureDict,
            std::shared_ptr<MeshRefiner> refinement_,
            markerMeshType modusFieldToMarker = markerMeshType::NonUniform,
            markerMeshType modusMarkerToField = markerMeshType::NonUniform
        );

        virtual void preSolveMovement()=0;
        virtual void preSolveMarkerMeshAdaption()=0;
        
        void solve() override;
        virtual void subTimestepStructureMovement()=0;
        virtual void subTimestepMarkerMeshAdaption()=0;
        
        // Adjoint methods
        void store() override;
        void setToTime(scalar time) override;
        
        // Sensitivity methods
        const volVectorField& getReferenceInOutField() const {return input_U;}
        const DynamicList<vector>& getMarkerCouplingForce() const {return markerCouplingForce;}
        const volVectorField& getVelocityField() const {return input_U;}
        
        void recomputeMarkerValues();
        
    protected:
        void interpolateFluidVelocityToMarkers();
        void computeCouplingForceOnMarkers();
        virtual void computeRodForceMoment();
        void interpolateFluidForceField();
        
        vector sumForces
        (
            std::function<bool(LagrangianMarker)> condition = [](LagrangianMarker marker){return true;}
        );
        
        vector sumMoments
        (
            std::function<bool(LagrangianMarker)> condition = [](LagrangianMarker marker){return true;}
        );
        
        bool printSummedRodForces = false;
        std::string recordRodForceFileName;
        std::unique_ptr<std::ofstream> recordRodForceFile;
        
        bool printSummedRodMoments = false;
        std::string recordRodMomentFileName;
        std::unique_ptr<std::ofstream> recordRodMomentFile;
        
        void moveMarkers();
        void meshMarkerAdaptation();
        
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

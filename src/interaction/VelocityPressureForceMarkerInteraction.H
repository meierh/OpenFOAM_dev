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
        
        void solve(bool finalIteration) override;
        virtual void subTimestepStructureMovement()=0;
        virtual void subTimestepMarkerMeshAdaption()=0;
        
        // Adjoint methods
        void store() override;
        void setToTime(scalar time) override;
        
        // Sensitivity methods
        const volVectorField& getReferenceInOutField() const {return input_U;}
        const DynamicList<vector>& getMarkerCouplingForce() const {return markerCouplingForce;}
        const volVectorField& getVelocityField() const {return input_U;}
        
        void recomputeMarkerValues() override;
        
        const DynamicList<vector>& getRodForce() const {return rodForce;}
        
    protected:
        void interpolateFluidVelocityToMarkers();
        void computeCouplingForceOnMarkers();
        virtual void computeRodForceMoment();
        void interpolateFluidForceField();

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
        
        class SumMarkerForceFile : public MarkerInfoFiles
        {
            public:
                SumMarkerForceFile
                (
                    const IOdictionary& structureDict
                ):
                MarkerInfoFiles
                (structureDict,"recordRodForce",true,{"Time","RodInd","Fx","Fy","Fz"})
                {}
            
                void writeSolution(const VelocityPressureForceInteraction& interaction);
        };        
        vector sumForces
        (
            std::function<bool(const LagrangianMarker&)> condition = [](const LagrangianMarker& marker){return true;}
        ) const ;
        SumMarkerForceFile sumMarkerForceFileObject;
        
        class SumMarkerMomentFile : public MarkerInfoFiles
        {
            public:
                SumMarkerMomentFile
                (
                    const IOdictionary& structureDict
                ):
                MarkerInfoFiles
                (structureDict,"recordRodMoment",true,{"Time","RodInd","Mx","My","Mz"})
                {}
            
                void writeSolution(const VelocityPressureForceInteraction& interaction);
        };
        vector sumMoments
        (
            std::function<bool(const LagrangianMarker&)> condition = [](const LagrangianMarker& marker){return true;}
        ) const;
        SumMarkerMomentFile sumMarkerMomentFileObject;
        
        class DetailedMarkerForceFile : public MarkerInfoFiles
        {
            public:
                DetailedMarkerForceFile
                (
                    const IOdictionary& structureDict
                ):
                MarkerInfoFiles
                (structureDict,"recordDetailedMarkerRodForce",true,{"Time","RodInd","Parameter","Angle","RadiusFrac","Px","Py","Pz","Nx","Ny","Nz","Fx","Fy","Fz"})
                {}
            
                void writeSolution(const VelocityPressureForceInteraction& interaction);
        };
        DetailedMarkerForceFile detailedMarkerForceFileObject;
};
}
#endif

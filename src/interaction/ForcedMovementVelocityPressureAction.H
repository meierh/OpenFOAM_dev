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
        
        ~ForcedMovementVelocityPressureAction();
               
        void preSolveMovement() override;
        void preSolveMarkerMeshAdaption() override;
        
        void subTimestepStructureMovement() override {};
        void subTimestepMarkerMeshAdaption() override {};
        
        void recomputeMarkerValues() override {NotImplemented;}
        
    protected:
        std::unique_ptr<List<List<vector>>> readDeformationDict();
        bool reconstructInterior(bool firstIteration,scalar time,bool finalIteration) override;
        vector getCellVelocity(label rodInd,scalar para,scalar angle,scalar radiusFrac) override;
        
    private:
        vector getVelocity(const LagrangianMarker* marker) override;
        std::vector<std::string> deformationNurbsCodes;
        void *rodMoveSo;
        std::vector<void(*)(unsigned int degree, const double* knots, unsigned int knotNbr, double* Px, double* Py, double* Pz, unsigned int PNbr, double time)> movementFunction;
        void constructMovementFunction();
        
        std::string rodMovementSignature = 
                "void deformationNurbs\
                (\
                    unsigned int degree, const double* knots, unsigned int knotNbr,\
                    double* Px, double* Py, double* Pz, unsigned int PNbr\
                )";
};

}
#endif

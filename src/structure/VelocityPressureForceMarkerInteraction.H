#ifndef VELOCITYPRESSUREFORCEMARKERINTERACTION_H
#define VELOCITYPRESSUREFORCEMARKERINTERACTION_H
#include "LineStructure.H"

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
            volVectorField& output_Uf
        );

        virtual void preSolveMovement()=0;
        
        void solve() override;
        
        void store() override;
        
        void setToTime(scalar time) override;
        
    protected:
        void interpolateFluidVelocityToMarkers();
        void computeCouplingForceOnMarkers();
        void computeRodForceMoment();
        void interpolateFluidForceField();
        
        vector sumForces
        (
            std::function<bool(LagrangianMarker)> condition = 
                                [](LagrangianMarker marker){return true;}
        );
        
        volVectorField& input_U;
        volVectorField& output_Uf;
        
        vector getVelocity(const LagrangianMarker* marker);
        
    private:
        DynamicList<vector> markerFluidVelocity;
        DynamicList<vector> makerCouplingForce;
        DynamicList<vector> rodForce;
        DynamicList<vector> rodMoment;
        
        std::map<scalar,std::tuple<DynamicList<vector>,DynamicList<vector>,DynamicList<vector>,DynamicList<vector>>> storage;
        
        void assignForceOnRod();
};
}
#endif

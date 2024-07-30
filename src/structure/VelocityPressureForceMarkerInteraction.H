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
            dynamicRefineFvMesh& mesh,
            LineStructure& structure,
            volVectorField& input_U,
            volVectorField& output_Uf
        );

        void interpolateFluidVelocityToMarkers();
        void computeCouplingForceOnMarkers();
        void computeRodForceMoment();
        void interpolateFluidForceField();
        
        vector sumForces
        (
            std::function<bool(LagrangianMarker)> condition = 
                                [](LagrangianMarker marker){return true;}
        );
        
    protected:
        volVectorField& input_U;
        volVectorField& output_Uf;
        
        vector getVelocity(const LagrangianMarker* marker);
        
    private:
        DynamicList<vector> markerFluidVelocity;
        DynamicList<vector> makerCouplingForce;
        DynamicList<vector> rodForce;
        DynamicList<vector> rodMoment;
        
        void assignForceOnRod();
};
}
#endif

#ifndef MARKERIMPLEMENTATION_H
#define MARKERIMPLEMENTATION_H
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
        
    private:
        volVectorField& input_U;
        volVectorField& output_Uf;
        
        DynamicList<vector> markerFluidVelocity;
        DynamicList<vector> makerCouplingForce;
        DynamicList<vector> rodForce;
        DynamicList<vector> rodMoment;
        
        void assignForceOnRod();
};

class TemperatureInteraction : public FieldMarkerStructureInteraction
{
    public:
        explicit TemperatureInteraction
        (
            dynamicRefineFvMesh& mesh,
            volScalarField& input_T,
            volScalarField& output_Tf
        );

        void interpolateTemperatureToMarkers();
        void computeCouplingHeatingOnMarkers();
        void computeRodHeating();
        void interpolateHeatingField();
        
        scalar sumHeating
        (
            std::function<bool(LagrangianMarker)> condition = 
                                [](LagrangianMarker marker){return true;}
        );
        
    private:
        volScalarField& input_T;
        volScalarField& output_Tf;
        
        DynamicList<scalar> markerFluidTemperature;
        DynamicList<scalar> makerCouplingHeating;
        DynamicList<scalar> rodHeating;
        
        void assignForceOnRod();
};

class TestInteraction : public FieldMarkerStructureInteraction
{
    public:
        explicit TestInteraction
        (
            dynamicRefineFvMesh& mesh,
            volScalarField& testField
        );
        
        void printSupportToField();
        
    private:
        volScalarField& testField;

};

}
#endif

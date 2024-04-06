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
            volVectorField& input_U,
            volScalarField& input_p,
            volVectorField& output_Uf
        );

        void interpolateFluidVelocityToMarkers();
        void computeCouplingForceOnMarkers();
        void interpolateFluidForceField();
        
    private:
        volVectorField& input_U;
        volScalarField& input_p;
        volVectorField& output_Uf;
        
        std::vector<vector> markerFluidVelocity;
        std::vector<vector> makerCouplingForce;
};
}
#endif

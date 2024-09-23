#ifndef FSIACTION_H
#define FSIACTION_H
#include "VelocityPressureForceMarkerInteraction.H"

namespace Foam
{

class FSIAction : public VelocityPressureForceInteraction
{
    public:
        explicit FSIAction
        (
            const fvMesh& mesh,
            LineStructure& structure,
            volVectorField& input_U,
            volVectorField& output_Uf,
            const IOdictionary& structureDict,
            std::shared_ptr<MeshRefiner> refinement_
        );
        
        void preSolveMovement() override {};
        
    protected:
        void computeRodForceMoment() override;
        
    private:
        vector getVelocity(const LagrangianMarker* marker) override;
        const IOdictionary& structureDict;
        std::pair<gsNurbs<scalar>,gsNurbs<scalar>> assignForceOnRod
        (
            const List<List<vector>>& globalRodForce,
            const List<List<vector>>& globalRodMoment,
            const List<List<label>>& globalRodNumber,
            const List<List<scalar>>& globalParameter
        );
};

}
#endif

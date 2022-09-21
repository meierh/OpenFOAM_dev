#ifndef FSIFLUIDBOUNDARYCONDITIONS_H
#define FSIFLUIDBOUNDARYCONDITIONS_H
#include "cutCellFvMesh.H"

namespace Foam
{
class FSIFluidBoundaryConditions
{
    public:
        FSIFluidBoundaryConditions
        (
            const dimensionedScalar& alpha,
            const volScalarField& T,
            const volScalarField& p,
            const volVectorField& U,
            const fvMesh& mesh,
            const dimensionedScalar nu,
            const word& IBpatchName = "cutCell"
        );
        FSIFluidBoundaryConditions(const FSIFluidBoundaryConditions&) = delete;
        
    private:
        const dimensionedScalar& alpha;
        const volScalarField& T;
        const volScalarField& p;
        const volVectorField& U;
        const fvMesh& mesh;
        label IBPatchID;
        const word& IBpatchName;
        const dimensionedScalar nu;
        const DynamicList<DynamicList<nurbsReference>>& meshPointNurbsReference
        
        void computeIBHeatFlux();
        void computeIBTemperature();
        void computeIBForce();
};
}

#endif

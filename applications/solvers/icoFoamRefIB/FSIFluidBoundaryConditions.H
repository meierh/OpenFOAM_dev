#ifndef FSIFLUIDBOUNDARYCONDITIONS_H
#define FSIFLUIDBOUNDARYCONDITIONS_H
#include "cutCellFvMesh.H"
#include "gismo.h"

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
            const cutCellFvMesh& mesh,
            const dimensionedScalar nu,
            const word& IBpatchName = "cutCell"
        );
        FSIFluidBoundaryConditions(const FSIFluidBoundaryConditions&) = delete;
        
    private:
        const dimensionedScalar& alpha;
        const volScalarField& T;
        const volScalarField& p;
        const volVectorField& U;
        const cutCellFvMesh& mesh;
        label IBPatchID;
        const word& IBpatchName;
        const dimensionedScalar nu;
        const DynamicList<DynamicList<cutCellFvMesh::nurbsReference>>& meshPointNurbsReference;
        const std::shared_ptr<std::vector<Nurbs1D>> Curves;
        std::unordered_map<label,label> boundaryPntToFace;
        List<std::multimap<scalar,label>> nurbsParameterToPnt;
        
        void computeIBHeatFlux();
        void computeIBTemperature();
        void computeIBForce();
        void assignBoundaryFacesToNurbsCurves();

};
}

#endif

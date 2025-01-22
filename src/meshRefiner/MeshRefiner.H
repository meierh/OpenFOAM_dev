#ifndef MESHREFINER_H
#define MESHREFINER_H
#include "fvMesh.H"
#include "LineStructure.H"
#include "fvMeshTopoChanger.H"
#include "fvMeshTopoChangersRefiner.H"

namespace Foam
{

class MeshRefiner
{
    public:
        explicit MeshRefiner
        (
            fvMesh& mesh,
            LineStructure& structure,
            volScalarField& doRefine,
            dictionary& dynamicMeshDict
        );
        virtual ~MeshRefiner() = default;
                
        /* preRefinedMesh: Was the mesh previously refined
         * Mesh is refined and markers are refined
         * Returns true if mesh was refined
         */
        bool refineMeshAndMarkers(bool preRefinedMesh=false);
        
        /* Refine only on fluid grounds
         * Returns true if mesh was refined
         */
        bool refineMeshOnFluid();
        
    protected:
        const scalar REFINE   = +1;
        const scalar MUSTKEEP = 0.25;
        const scalar DONTCARE = -0.25;
        const scalar UNREFINE = -1;

        const scalar REFINE_LIM   = +0.5;
        const scalar MUSTKEEP_LIM =  0;
        const scalar DONTCARE_LIM = -0.5;

        fvMesh& mesh;
        LineStructure& structure;
        volScalarField& doRefine;
        dictionary& dynamicMeshDict;

        volScalarField fieldRefineDemands;
        virtual void fieldRefinement();

        volScalarField markerRefineDemands;
        void markerRefinement(scalar defaultValue);
        
        bool applyMeshAdaption();
        
        const fvMeshTopoChangers::refiner* refinerPtr;

    private:
        scalar refinementDemandMerge(scalar fieldDemand, scalar markerDemand);
        scalar markerCharLengthToCellSizeFactor = 0.1;
};

}
#endif

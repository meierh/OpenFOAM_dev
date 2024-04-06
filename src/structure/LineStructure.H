#ifndef LINESTRUCTURE_H
#define LINESTRUCTURE_H
#include "Structure.H"
//#include "rodTools.h"

namespace Foam
{

class LagrangianMarker
{
    public:
        explicit LagrangianMarker
        (
            const scalar markerParameter,
            const vector markerPosition,
            const ActiveRodMesh::rodCosserat* baseRod
        ):
        markerParameter(markerParameter),
        markerPosition(markerPosition),
        baseRod(baseRod)
        {}
        
        vector getMarkerVelocity();

        const scalar markerParameter;
        const vector markerPosition;
        const ActiveRodMesh::rodCosserat* baseRod;
        scalar markerVolume;
        std::unordered_set<label> supportCells;
        label markerCell;
        std::string to_string();
        
};

class FieldMarkerStructureInteraction
{
    public:
        explicit FieldMarkerStructureInteraction
        (
            dynamicRefineFvMesh& mesh
        );
        
        ActiveRodMesh::rodMesh* rodMesh;
        std::vector<LagrangianMarker*> markers;
        dynamicRefineFvMesh& mesh;
        
    protected:
        scalar phiFunction(scalar r);
        scalar deltaDirac(vector X, vector x, scalar h);
        template<typename T>
        void fieldToMarker
        (
            GeometricField<T,fvPatchField,volMesh>& fieldData,
            std::vector<T>& markerData
        )
        {
            Info<<"Transfer field to markers"<<Foam::endl;
            const cellList& cells = mesh.cells();
            const faceList& faces = mesh.faces();
            const pointField& points = mesh.points();
            
            markerData.resize(markers.size());
            for(label index=0; index<markerData.size(); index++)
            {
                T value = Foam::zero();
                const LagrangianMarker& oneMarker = *(markers[index]);
                const vector& oneMarkerPos = oneMarker.markerPosition;
                for(auto cellIter=oneMarker.supportCells.begin();
                    cellIter!=oneMarker.supportCells.end();
                    cellIter++)
                {
                    label cellInd = *cellIter;
                    const cell& oneCell = cells[cellInd];
                    scalar cellVol = oneCell.mag(points,faces);
                    point cellCentre = oneCell.centre(points,faces);
                    scalar h = std::cbrt(cellVol);
                    scalar factor = deltaDirac(cellCentre,oneMarkerPos,h)*cellVol;
                    value += fieldData[cellInd]*factor;
                }
                markerData[index] = value;
            }
        }
        template<typename T>
        void markerToField
        (
            std::vector<T>& markerData,
            GeometricField<T,fvPatchField,volMesh>& fieldData
        )
        {
            const cellList& cells = mesh.cells();
            const faceList& faces = mesh.faces();
            const pointField& points = mesh.points();
            
            if(markerData.size()!=markers.size())
                FatalErrorInFunction<<"Marker Data size mismatch"<< exit(FatalError);
            
            List<DynamicList<label>> cellMarkers(cells.size());
            for(label markerInd=0; markerInd<markers.size(); markerInd++)
            {
                const LagrangianMarker& oneMarker = *(markers[markerInd]);
                for(auto cellIter=oneMarker.supportCells.begin();
                    cellIter!=oneMarker.supportCells.end();
                    cellIter++)
                {
                    label cellInd = *cellIter;
                    cellMarkers[cellInd].append(markerInd);
                }
            }
            
            for(label cellInd=0; cellInd<fieldData.size(); cellInd++)
            {
                const cell& oneCell = cells[cellInd];
                point cellCentre = oneCell.centre(points,faces);
                scalar cellVol = oneCell.mag(points,faces);
                scalar h = std::cbrt(cellVol);
                if(cellMarkers[cellInd].size()>0)
                {
                    T value = Foam::zero();
                    for(label markerInd : cellMarkers[cellInd])
                    {
                        const LagrangianMarker& oneMarker = *(markers[markerInd]);
                        const vector& oneMarkerPos = oneMarker.markerPosition;
                        scalar markerVolume = oneMarker.markerVolume;
                        scalar factor = deltaDirac(cellCentre,oneMarkerPos,h)*markerVolume;
                        value += markerData[markerInd]*factor;
                    }
                    fieldData[cellInd] = value;
                }
                else
                {
                    fieldData[cellInd] = Foam::zero();
                }
            }
        }
};

class LineStructure : public Structure
{
    public:
        explicit LineStructure
        (
            dynamicRefineFvMesh& mesh,
            const dimensionedScalar& alpha,
            const volScalarField& T,
            const volScalarField& p,
            const volVectorField& U,
            const dimensionedScalar nu,
            const List<scalar> crossSecArea
        );
        LineStructure(const LineStructure&) = delete;
        void transferMarkers(FieldMarkerStructureInteraction& connector);

    protected:
        void getSupportDomain
        (
            label cellInd,
            std::unordered_set<label>& neighbourhoodCells
        );
        
        std::pair<scalar,scalar> minMaxSpan (const cell& thisCell);
        
        void cellDistances
        (
            const cell& thisCell,
            std::vector<scalar>&
        );
        
        scalar supportDomainMinSize(label cellInd);
        scalar supportDomainMinSize(const std::unordered_set<label>& supportDomainCells);
        
        static scalar evaluateRodArcLen
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parA,
            scalar parB
        );
        
        static vector evaluateRodPos
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parameter
        );
        
        /*
        vector evaluateRodDeriv
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parameter
        );
        */
        
        bool doSubdivision
        (
            const LagrangianMarker& smallerSide,
            const LagrangianMarker& largerSide
        );
        
        LagrangianMarker createLagrangianMarker
        (
            scalar markerParameter,
            const ActiveRodMesh::rodCosserat* oneRod
        );
        
        static scalar distance(const LagrangianMarker& A, const LagrangianMarker& B);
        
        const label minRefinement = 5;
        
    private:
        const List<scalar> crossSecArea;
        
        std::vector<std::unique_ptr<std::vector<LagrangianMarker>>> rodMarkers;
        
        std::unique_ptr<std::vector<LagrangianMarker>> constructMarkerSet
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar crossSecArea
        );
        
};
}
#endif

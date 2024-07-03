#ifndef LINESTRUCTURE_H
#define LINESTRUCTURE_H
#include <iterator>
#include "Structure.H"
#include "LinearSystem.H"

namespace Foam
{

class LagrangianMarker
{
    public:
        explicit LagrangianMarker
        (
            const Structure& structure,
            const dynamicRefineFvMesh& mesh,
            const label rodNumber,
            const ActiveRodMesh::rodCosserat* baseRod,
            const scalar markerParameter
        );
    protected:
        explicit LagrangianMarker
        (    
            const Structure& structure,
            const dynamicRefineFvMesh& mesh,
            const label rodNumber,
            const ActiveRodMesh::rodCosserat* baseRod
        );
        
    public:
        vector getMarkerVelocity();
        scalar getMarkerTemperature();
        scalar getMarkerCellVolume();
        std::pair<scalar,scalar> getMarkerCellSpacing();
        scalar getMarkerCellMinSpacing();
        std::string to_string() const;
        void evaluateMarker();
        
        inline label getRodNumber() const {return rodNumber;}
        inline const ActiveRodMesh::rodCosserat* getBaseRod() const {return baseRod;}
        inline scalar getMarkerParameter() const {return markerParameter;}
        inline vector getMarkerPosition() const {return markerPosition;}
        inline label getMarkerCell() const {return markerCell;}
        inline void setMarkerVolume(scalar volume) {markerVolume=volume;}
        inline scalar getMarkerVolume() const {return markerVolume;}
        inline void setMarkerWeight(scalar weight) {markerWeight=weight;}
        inline scalar getMarkerWeight() const {return markerWeight;}
        inline const DynamicList<std::tuple<bool,label,label>>& getSupportCells() const {return supportCells;}
        inline vector getDilation() const {return dilation;}
        inline const std::array<scalar,10>& getCorrParaB() const {return b;}
        
        static scalar deltaDirac(vector X, vector x, scalar h);
        static scalar deltaDirac(vector X, vector x, vector h);
        static scalar correctedDeltaDirac(vector X, vector x, scalar h, const std::array<scalar,10>& b);
        static scalar correctedDeltaDirac(vector X, vector x, vector h, const std::array<scalar,10>& b);
        static scalar phiFunction(scalar r);
        
        //h=getDilation() b=b;
        scalar correctedDeltaDirac(vector X, vector x) const;
        scalar deltaDirac(vector X, vector x) const;
        
        void compNonUniformCorrWeights() {computeCorrectionWeights();}
        
        template<typename VALTYPE>
        VALTYPE convolute
        (
            std::function<scalar(vector,vector)> weightFunction,
            std::function<VALTYPE(std::tuple<bool,label,label>)> valueFunction
        ) const
        {
            vector X = getMarkerPosition();
            const DynamicList<std::tuple<bool,label,label>>& supportCells = getSupportCells();
            
            VALTYPE convValue = Foam::zero();
            for(label i=0; i<supportCells.size(); i++)
            {
                const std::tuple<bool,label,label>& suppCellData = supportCells[i];
                vector cellCentre;
                scalar cellVolume;
                getCellData(suppCellData,cellCentre,cellVolume);
                vector x = cellCentre;

                scalar w = weightFunction(X,x);
                VALTYPE f = valueFunction(suppCellData);
                
                convValue += w*f*cellVolume;
            }
            return convValue;
        }
        void getCellData(const std::tuple<bool,label,label>& cell, vector& cellCentre, scalar& volume) const;
        
    protected:        
        const Structure& structure;
        const dynamicRefineFvMesh& mesh;
        label rodNumber;
        const ActiveRodMesh::rodCosserat* baseRod;
        scalar markerParameter;
        vector markerPosition;
        label markerCell;
        scalar markerVolume;
        scalar markerWeight = 1;
        
        // std::tuple<own processor,processor ID,cell ID>
        DynamicList<std::tuple<bool,label,label>> supportCells;
        Vector<bool> existingDims;
        vector h_plus;
        vector h_minus;
        vector dilation;
        
        void computeSupport(label iterations=1);
        void computeSupportDimensions();
        void minMaxSupportWidth();
        void dilationFactors();
        
        std::array<scalar,10> b;
        
        scalar computeMoment(vector indices, std::function<scalar(vector,vector)> deltaFunction) const;
        scalar computeCorrectedMoment(vector indices) const;
        scalar computeMoment(vector indices) const;
        std::unique_ptr<gismo::gsMatrix<scalar>> computeMomentMatrix3D() const ;
        std::unique_ptr<gismo::gsMatrix<scalar>> computeMomentMatrix2D(std::array<label,2> dim) const;
        std::unique_ptr<gismo::gsMatrix<scalar>> computeMomentMatrix1D(label dim) const;
        std::unique_ptr<std::array<scalar,10>> rescalingDiagonal3D() const;
        std::unique_ptr<std::array<scalar,6>> rescalingDiagonal2D(label normalDim) const;
        std::unique_ptr<std::array<scalar,3>> rescalingDiagonal1D(label dim) const;
        void computeCorrectionWeights();
        
        enum WeightSystemSolve {Raw,GeoRescaled,DiagRescaled};
        template<uint dim>
        void solveWeights
        (
            std::unique_ptr<std::array<scalar,dim>> Hptr,
            std::unique_ptr<gismo::gsMatrix<scalar>> Mptr,
            gismo::gsMatrix<scalar>& x,
            WeightSystemSolve solutionStrategy = WeightSystemSolve::GeoRescaled
        )
        {
            std::array<scalar,dim>& H = *Hptr;
            gismo::gsMatrix<scalar>& M = *Mptr;
            x = gismo::gsMatrix<scalar>(dim,1);
            gismo::gsMatrix<scalar> e(dim,1);
            for(uint i=0; i<dim; i++)
            {
                e(i,0) = 0;
                x(i,0) = 0;
            }
            e(0,0) = 1;
            
            gismo::gsMatrix<scalar> P(dim,dim);
            for(uint r=0;r<dim;r++)
                for(uint c=0;c<dim;c++)
                    P(r,c)=0;
            if(solutionStrategy==WeightSystemSolve::Raw)
            {
                for(uint d=0; d<dim; d++)
                    P(d,d) = 1;
            }
            else if(solutionStrategy==WeightSystemSolve::GeoRescaled)
            {
                for(uint d=0; d<dim; d++)
                    P(d,d) = H[d];
            }
            else
            {
                for(uint d=0; d<dim; d++)
                    P(d,d) = 1.0/M(d,d);
            }
            gismo::gsMatrix<scalar> A = P*M;

            linearSolve(A,x,e);
        }

        std::unique_ptr<gismo::gsMatrix<scalar>> computeCorrectedMomentMatrix() const;
};

class LineStructure;

class FieldMarkerStructureInteraction
{
    public:
        enum markerMeshType {Uniform, NonUniform};
        explicit FieldMarkerStructureInteraction
        (
            dynamicRefineFvMesh& mesh,
            LineStructure& structure,
            markerMeshType modusFieldToMarker = markerMeshType::Uniform,
            markerMeshType modusMarkerToField = markerMeshType::Uniform
        );
        ActiveRodMesh::rodMesh* rodMesh;
        std::vector<LagrangianMarker*> markers;
        
        dynamicRefineFvMesh& mesh;
        LineStructure& structure;
        
    protected:

        //Map field from fielData to markerData according to markers
        template<typename T> void fieldToMarker
        (
            const GeometricField<T,fvPatchField,volMesh>& fieldData,
            DynamicList<T>& markerData
        );
                
        template<typename T> void markerToField
        (
            const DynamicList<T>& markerData,
            GeometricField<T,fvPatchField,volMesh>& fieldData
        );
        
        void scatterNurbs(std::pair<gsNurbs<scalar>,label> in, std::pair<gsNurbs<scalar>,label>& out);
        
        template<typename T> void gatherMarkerData
        (
            std::tuple<List<scalar>,List<List<T>>,label> in,
            std::tuple<List<scalar>,List<List<T>>,label>& out
        );
        
        const scalar h;
        const markerMeshType modusFieldToMarker;
        const markerMeshType modusMarkerToField;
};

class LineStructure : public Structure
{
    public:
        enum markerMeshType {Uniform=0, NonUniform=1};
        explicit LineStructure
        (
            dynamicRefineFvMesh& mesh,
            const List<scalar> crossSecArea,            
            markerMeshType modusFieldToMarker = markerMeshType::NonUniform,
            markerMeshType modusMarkerToField = markerMeshType::NonUniform
        );
        LineStructure(const LineStructure&) = delete;
        void connect(FieldMarkerStructureInteraction& connector);

        void createMarkers();

        void reInitializeMarkers(bool keepMarkers=false,bool keepSeedPoints=false);

        static vector evaluateRodPos
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parameter
        );

        void test();

    protected:
        
        void createSpacingPoints();
        void createSpacedPointsOnRod
        (
            label rodNumber,
            scalar spacing
        );

        void createMarkersFromSpacedPoints();
        void createMarkersFromSpacedPointsOnRod
        (
            label rodNumber
        );

        void refineMarkers
        (
            std::pair<bool,scalar> forcedSpacing = {false,-1}
        );
        void refineMarkersOnRod
        (
            label rodNumber,
            std::pair<bool,scalar> forcedSpacing = {false,-1}
        );

        void setMarkerVolume();
        void setMarkerVolumeOnRod
        (
            label rodNumber
        );

        void evaluateMarkerMeshRelation();
        void evaluateMarkerMeshRelation(std::list<LagrangianMarker>& markerList);

        void reduceMarkers();

        template<typename MARKER>
        class MarkerReference
        {
            public:
                MarkerReference
                (
                    typename std::list<MARKER>::iterator iterator,
                    typename std::list<MARKER>* listPointer
                ):
                iterator(iterator),listPointer(listPointer),markedToDelete(false)
                {}
                
                const MARKER& getMarker() const {return *iterator;}
                MARKER& getMarker() {return *iterator;}
                void deleteMarker() {markedToDelete=true;}
                
                ~MarkerReference()
                {
                    if(markedToDelete)
                    {
                        listPointer->erase(iterator);
                    }
                }
            
            private:
                typename std::list<MARKER>::iterator iterator;
                typename std::list<MARKER>* listPointer;
                bool markedToDelete;
        };

        template<class MARKER>
        void reduceMarkers
        (
            std::vector<MarkerReference<MARKER>>& allMarkers
        )
        {
            using SingleMarker = MarkerReference<MARKER>;
            
            std::unordered_map<label,std::vector<SingleMarker*>> cellToMarker;
            for(SingleMarker& markerData : allMarkers)
            {
                const MARKER& marker = markerData.getMarker();
                label cellInd = marker.getMarkerCell();
                cellToMarker[cellInd].push_back(&markerData);
            }
            oneCenteredMarker(allMarkers,cellToMarker);
        }

        template<class MARKER>
        void maxOneMarkerPerOctant
        (
            std::vector<MarkerReference<MARKER>>& allMarkers,
            std::unordered_map<label,std::vector<MarkerReference<MARKER>*>>& cellToMarker
        );

        template<class MARKER>
        void maxOneMarkerPerCell
        (
            std::vector<MarkerReference<MARKER>>& allMarkers,
            std::unordered_map<label,std::vector<MarkerReference<MARKER>*>>& cellToMarker
        );

        template<class MARKER>
        void oneCenteredMarker
        (
            std::vector<MarkerReference<MARKER>>& allMarkers,
            std::unordered_map<label,std::vector<MarkerReference<MARKER>*>>& cellToMarker
        );

        void collectMarkers();
        
        void collectHaloMarkers();
        
        void exchangeHaloMarkersData();
                
        void computeMarkerCellWeights();

        using LinearSystem = std::tuple<CSR_Matrix_par,Vector_par>;
        std::unique_ptr<LinearSystem> computeMarkerEpsilonMatrix();
        void computeMarkerWeights();
                
        static scalar evaluateRodArcLen
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parA,
            scalar parB
        );
                    
        template<typename T>
        static T integrateRodwise
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            std::function<T(scalar)> function
        )
        {
            return integrateRodwise<T>(oneRod,oneRod->m_Curve.knots().first(),oneRod->m_Curve.knots().last(),function);
        }

        template<typename T>
        static T integrateRodwise
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parA,
            scalar parB,
            std::function<T(scalar)> function
        )
        {
            if(parA>parB)
                FatalErrorInFunction<<"Integration from higher to lower value"<<exit(FatalError);
            
            std::set<scalar> knotSet;
            for(scalar knot : oneRod->m_Curve.knots())
                knotSet.insert(knot);
            
            label degree = oneRod->m_Curve.knots().degree();
            auto start = knotSet.insert(parA);
            auto end = knotSet.insert(parB);
            
            label numberOfAbcissa = (degree+1)*(degree+1);
            
            T totalValue = Foam::zero();
            
            std::set<scalar>::iterator currNode = start.first;
            std::set<scalar>::iterator nextNode = ++(start.first);
            for(; currNode!=knotSet.end() && nextNode!=knotSet.end() && currNode!=end.first ; )
            {
                scalar startPar = *currNode;
                scalar endPar = *nextNode;
                scalar dist = endPar-startPar;
                scalar step = dist/numberOfAbcissa;
                scalar initialStep = startPar+step/2;
                T summedValue = Foam::zero();
                for(label i=0;i<numberOfAbcissa;i++)
                {
                    summedValue += function(initialStep)*step;
                    initialStep+=step;
                }
                totalValue+=summedValue;
                
                currNode=nextNode;
                nextNode++;
            }
            return totalValue;
        }
        
        static scalar distance
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parA,
            scalar parB
        );
        
        static scalar distance(const LagrangianMarker& A, const LagrangianMarker& B);
        
        std::vector<FieldMarkerStructureInteraction*> connectedInteractions;
        //const label minRefinement = 1;
        markerMeshType modusFieldToMarker;
        markerMeshType modusMarkerToField;
        
    private:
        const List<scalar> crossSecArea;
        // [rodNumber] -> ...
    protected:
        /********************Marker points**************************/
        std::vector<std::unique_ptr<std::vector<scalar>>> initialRodPoints;
        std::vector<std::unique_ptr<std::list<LagrangianMarker>>> rodMarkersList;
        std::vector<LagrangianMarker*> collectedMarkers;
        
        /******************Halo markers data************************/
        class GlobalHaloMarkers
        {
            public:
                GlobalHaloMarkers
                (
                    std::vector<std::vector<std::pair<LagrangianMarker*,label>>>& selfhaloCellsRodMarkersList
                );
                GlobalHaloMarkers(){}
                
                void appendMarkerData
                (
                    label haloCellInd,
                    std::tuple<vector,scalar,label,vector,DynamicList<Pair<label>>,DynamicList<vector>,DynamicList<scalar>,List<scalar>> data
                );
                
                std::tuple<vector,scalar,label,vector,DynamicList<Pair<label>>,DynamicList<vector>,DynamicList<scalar>,List<scalar>> getMarkerData
                (
                    label process,
                    label haloCellInd,
                    label cellMarkerInd
                ) const;
                
                inline label size_processes() const;
                inline label size_haloCells(label process) const;
                inline label size_cellMarkers(label process, label haloCellInd) const;
                
                void communicate();
                
                template<typename T>
                std::unique_ptr<List<List<DynamicList<T>>>> broadcasHaloMarkerField
                (
                    const List<T>& collectedMarkerField
                );
            
            private:
                // [haloCells] -> [marker of cell] -> ptr,index
                std::vector<std::vector<std::pair<LagrangianMarker*,label>>>* selfhaloCellsRodMarkersList;
                
                label nProcs;
                List<label> procHaloCellsSize;
                List<List<label>> procHaloCellMarkerSize;
                
                // [nProc] -> [haloCells] -> [markers] -> position
                List<List<DynamicList<vector>>> globalHaloCellsMarkerPos;
                // [nProc] -> [haloCells] -> [markers] -> volume
                List<List<DynamicList<scalar>>> globalHaloCellsMarkerVolume;
                // [nProc] -> [haloCells] -> [markers] -> localIndex
                List<List<DynamicList<label>>> globalHaloCellsLocalIndex;
                // [nProc] -> [haloCells] -> [markers] -> dilation
                List<List<DynamicList<vector>>> globalHaloCellsMarkerDilation;
                // [nProc] -> [haloCells] -> [markers] -> supportCellIndices{proc,cellInd}[]
                List<List<DynamicList<DynamicList<Pair<label>>>>> globalHaloCellsMarkerSupportCellIndices;
                // [nProc] -> [haloCells] -> [markers] -> supportCellsCentre[]
                List<List<DynamicList<DynamicList<vector>>>> globalHaloCellsMarkerSupportCellCentres;
                // [nProc] -> [haloCells] -> [markers] -> supportCellsVolume[]
                List<List<DynamicList<DynamicList<scalar>>>> globalHaloCellsMarkerSupportCellVolume;
                // [nProc] -> [haloCells] -> [markers] -> b[]
                List<List<DynamicList<List<scalar>>>> globalHaloCellsMarkerb;
                
                void check();
                void check(label process, label haloCellInd, label cellMarkerInd) const;
                bool broadcasted;
        };
        // [haloCells] -> [marker of cell] -> ptr        
        std::vector<std::vector<std::pair<LagrangianMarker*,label>>> haloCellsRodMarkersList;
        std::unique_ptr<List<List<DynamicList<scalar>>>> haloMarkerWeights;
        GlobalHaloMarkers globHaloMarkers;
    public:
        inline std::unique_ptr<List<List<DynamicList<scalar>>>>& get_haloMarkerWeights(){return haloMarkerWeights;}
        inline const GlobalHaloMarkers& get_globalHaloMarkers(){return globHaloMarkers;}
        
        template<typename T>
        std::unique_ptr<List<List<DynamicList<scalar>>>> broadcastHaloMarkerFields
        (
            const List<T>& collectedMarkerField
        );
        
        scalar initialSpacingFromMesh();
        void check();
};
#include "LineStructure.T.H"
}
#endif

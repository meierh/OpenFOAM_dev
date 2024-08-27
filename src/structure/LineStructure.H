#ifndef LINESTRUCTURE_H
#define LINESTRUCTURE_H
#include <iterator>
#include "Structure.H"
#include "LinearSystem.H"

namespace Foam
{

enum SystemSolve {Raw,RowAequilibration,ColAequilibration,GeoRescaled,Jacobi};   
    
class LagrangianMarker
{
    public:
        enum SupportType {Direct,Full};
        virtual ~LagrangianMarker()=default;
        
        explicit LagrangianMarker
        (
            const Structure& structure,
            const fvMesh& mesh,
            const label rodNumber,
            const ActiveRodMesh::rodCosserat* baseRod,
            const scalar markerParameter
        );
    protected:
        explicit LagrangianMarker
        (    
            const Structure& structure,
            const fvMesh& mesh,
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
        void total_print() const;
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
        inline scalar getMarkerCharacLen() const {return markerCharLen;}
        virtual inline scalar getMarkerAngle() const {return 0;}
        virtual inline scalar getMarkerRadiusFrac() const {return 0;}
        inline const DynamicList<Pair<label>>& getSupportCells(SupportType type=SupportType::Direct) const
        {
            if(type==SupportType::Direct)
                return directSupport;
            else
                return fullSupport;
        }
        inline vector getDilation() const {return dilation;}
        inline const FixedList<scalar,10>& getCorrParaB() const {return b;}
        
        static scalar deltaDirac(vector X, vector x, scalar h);
        static scalar deltaDirac(vector X, vector x, vector h);
        static scalar correctedDeltaDirac(vector X, vector x, scalar h, const FixedList<scalar,10>& b);
        static scalar correctedDeltaDirac(vector X, vector x, vector h, const FixedList<scalar,10>& b);
        static scalar phiFunction(scalar r);
        
        static scalar ddeltaDirac_dX(vector X, vector x, vector h);
        static scalar dcorrectedDeltaDirac_dX(vector X, vector x, vector h, const FixedList<scalar,10>& b);
        
        //h=getDilation() b=b;
        scalar correctedDeltaDirac(vector X, vector x) const;
        scalar deltaDirac(vector X, vector x) const;
        
        vector dcorrectedDeltaDirac_dX(vector X, vector x) const;
        vector ddeltaDirac_dX(vector X, vector x) const;
        
        void compNonUniformCorrWeights() {computeCorrectionWeights();}
        
        template<typename VALTYPE>
        VALTYPE convolute
        (
            std::function<scalar(vector,vector)> weightFunction,
            std::function<VALTYPE(Pair<label>)> valueFunction
        ) const;
        
        template<typename VALTYPE>
        VALTYPE convolute
        (
            std::function<scalar(vector,vector)> weightFunction,
            std::function<VALTYPE(Pair<label>)> valueFunction,
            const List<Pair<label>>& support
        ) const;
        
        void getCellData(const Pair<label>& cell, vector& cellCentre, scalar& volume) const;
        
    protected:        
        const Structure& structure;
        const fvMesh& mesh;
        label rodNumber;
        const ActiveRodMesh::rodCosserat* baseRod;
        scalar markerParameter;
        vector markerPosition;
        label markerCell;
        scalar markerVolume;
        scalar markerWeight = 1;
        scalar markerCharLen;
        
        // Pair<processor,cell ID>
        const label supportWidth = 2;
        DynamicList<Pair<label>> directSupport;
        DynamicList<Pair<label>> fullSupport;
        SupportType momentsSupportType = SupportType::Full;
        Vector<bool> existingDims;
        void computeSupport(label iterations=2);
        
        vector h_plus;
        vector h_minus;
        Pair<vector> minMaxNeighbourWidth(const List<Pair<label>>& support) const;
        
        vector dilation;
        vector dilationFactors(Pair<vector> h) const;
        
        void checkDirectSupport() const;
        void reduceSupport();
        
        FixedList<scalar,10> b;
        
        scalar computeMoment(vector indices, std::function<scalar(vector,vector)> deltaFunction) const;
        scalar computeCorrectedMoment(vector indices) const;
        scalar computeMoment(vector indices) const;
        std::unique_ptr<gismo::gsMatrix<scalar>> computeMomentMatrix3D() const ;
        std::unique_ptr<gismo::gsMatrix<scalar>> computeMomentMatrix2D(std::array<label,2> dim) const;
        std::unique_ptr<gismo::gsMatrix<scalar>> computeMomentMatrix1D(label dim) const;
        std::unique_ptr<gismo::gsMatrix<scalar>> rescalingDiagonal3D() const;
        std::unique_ptr<gismo::gsMatrix<scalar>> rescalingDiagonal2D(label normalDim) const;
        std::unique_ptr<gismo::gsMatrix<scalar>> rescalingDiagonal1D(label dim) const;
        std::unique_ptr<Pair<gismo::gsMatrix<scalar>>> computeMomentsMatrix(Vector<bool> dims) const;
        std::unique_ptr<Pair<gismo::gsMatrix<scalar>>> computeMomentsMatrix
        (
            Vector<bool> dims,
            SystemSolve solutionStrategy
        ) const;
        Vector<bool> analyseMomentsMatrix(Vector<bool> dimensions, const gismo::gsMatrix<scalar>& momentsMatrix) const;
        bool checkSolvability
        (
            std::unique_ptr<Pair<gismo::gsMatrix<scalar>>>& system,
            Vector<bool> dimensions
        );
        void searchValidConvolutionSetup(std::unique_ptr<Pair<gismo::gsMatrix<scalar>>>& system);
        
        void computeCorrectionWeights();
        
        SystemSolve solutionStrategy = SystemSolve::GeoRescaled;
        scalar conditionThreshold = 1e5;
        std::unique_ptr<gismo::gsMatrix<scalar>> computeCorrectedMomentMatrix() const;
};
#include "LagrangianMarker.T.H"

class LineStructureStatus
{
    public:
        class Node
        {
            public:
                Node(bool value, std::string name):value(value),name(name){}
                Node(std::vector<Node*> dependencies, bool value, std::string name):dependencies(dependencies),value(value),name(name){}
                std::vector<Node*> dependencies;
                bool value;
                std::vector<Node*> following;
                std::string name;
        };
        LineStructureStatus();
        void execValid(const Node& method);
        void executed(Node& method);

        Node initialPoints =        Node(false,"initialPoints");
        Node markers =              Node({&initialPoints},false,"markers");
        Node markersRefined =       Node({&markers},false,"markersRefined");
        Node markersVolume =        Node({&markers,&markersRefined},false,"markersVolume");
        Node markersReduction =     Node({&markers,&markersVolume},false,"markersReduction");
        Node markerMesh =           Node({&markers},false,"markerMesh");
        Node markersCollected =     Node({&markers},false,"markersCollected");
        Node markersCellWeight =    Node({&markersCollected,&markerMesh},false,"markersCellWeight");
        Node markersHaloCollect =   Node({&markersCollected},false,"markersHaloCollect");
        Node markersHaloExchange =  Node({&markersHaloCollect},false,"markersHaloExchange");
        Node markersWeight =        Node({&markersCellWeight,&markersHaloExchange},false,"markersWeight");
        Node markersWeightExchange= Node({&markersWeight},false,"markersWeightExchange");
};

class LineStructure;

class FieldMarkerStructureInteraction
{
    public:
        enum markerMeshType {Uniform, NonUniform};
        explicit FieldMarkerStructureInteraction
        (
            fvMesh& mesh,
            LineStructure& structure,
            markerMeshType modusFieldToMarker = markerMeshType::NonUniform,
            markerMeshType modusMarkerToField = markerMeshType::NonUniform
        );
        ActiveRodMesh::rodMesh* rodMesh;
        //std::vector<LagrangianMarker*> markers;
        
        fvMesh& mesh;
        LineStructure& structure;
        
    protected:

        //Map field from fielData to markerData according to markers
        template<typename T> void fieldToMarker
        (
            const GeometricField<T,fvPatchField,volMesh>& fieldData,
            DynamicList<T>& markerData
        );
        
        template<typename T> void integrateFieldOverMarker
        (
            const GeometricField<T,fvPatchField,volMesh>& fieldData,
            DynamicList<T>& markerData,
            std::function<T(const LagrangianMarker*, std::function<T(Pair<label>)>)> integration            
        );
                
        template<typename T> void markerToField
        (
            const DynamicList<T>& markerData,
            GeometricField<T,fvPatchField,volMesh>& fieldData
        );
                
        template<typename T> void integrateMarkersOverField
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

class NurbsCoeffReference
{
    public:
        NurbsCoeffReference(label rodNumber, label coeffNumber, label dimension);
        label rodNumber;
        label coeffNumber;
        label dimension;
};

struct NurbsCoeffHash
{
    label operator()(const NurbsCoeffReference &r) const
    {
        return std::hash<label>{}(r.rodNumber)^std::hash<label>{}(r.coeffNumber)^std::hash<label>{}(r.dimension);
    }
};

class CrossSectionCoeffReference
{
    public:
        CrossSectionCoeffReference(label rodNumber, label fourierCoeffNumber, label coeffNumber);
        label rodNumber;
        label fourierCoeffNumber;
        label coeffNumber;
};

struct CrossSecCoeffHash
{
    label operator()(const CrossSectionCoeffReference &r) const
    {
        return std::hash<label>{}(r.rodNumber)^std::hash<label>{}(r.fourierCoeffNumber)^std::hash<label>{}(r.coeffNumber);
    }
};

class Parameter
{
    public:
        enum Type {Rod,CrossSection,None};   
        Parameter();
        
        Parameter(NurbsCoeffReference coeff);
        Parameter(const std::vector<NurbsCoeffReference>& coeffs);
        void addCoeff(NurbsCoeffReference coeffRef);
        
        Parameter(CrossSectionCoeffReference coeff);
        Parameter(const std::vector<CrossSectionCoeffReference>& coeffs);
        void addCoeff(CrossSectionCoeffReference ref);
        
        label getDimension() const {return dimension;}
        bool isValid() const {return valid;}
        Type getType() const {return parameterType;}
        
        const std::vector<NurbsCoeffReference>& getNurbsCoeffs() const {return nurbsCoeffs;}
        const std::vector<CrossSectionCoeffReference>& getCrossSecCoeffs() const {return crossSecCoeffs;}
        
    private:
        Type parameterType;
        
        std::vector<NurbsCoeffReference> nurbsCoeffs;
        label dimension;
        
        std::vector<CrossSectionCoeffReference> crossSecCoeffs;
        
        bool valid=false;
};

class LineStructure : public Structure
{
    public:
        enum markerMeshType {Uniform=0, NonUniform=1};
        explicit LineStructure
        (
            fvMesh& mesh,
            const List<scalar> crossSecArea,            
            markerMeshType modusFieldToMarker = markerMeshType::NonUniform,
            markerMeshType modusMarkerToField = markerMeshType::NonUniform
        );
        explicit LineStructure
        (
            fvMesh& mesh,
            markerMeshType modusFieldToMarker = markerMeshType::NonUniform,
            markerMeshType modusMarkerToField = markerMeshType::NonUniform
        );
        explicit LineStructure
        (
            fvMesh& mesh,
            const scalar crossSecArea,
            markerMeshType modusFieldToMarker = markerMeshType::NonUniform,
            markerMeshType modusMarkerToField = markerMeshType::NonUniform
        );
        LineStructure(const LineStructure&) = delete;
        virtual ~LineStructure() = default;
        
        void connect(FieldMarkerStructureInteraction& connector);

        void createMarkers();

        void reInitializeMarkers(bool keepMarkers=false,bool keepSeedPoints=false);

        static vector evaluateRodPos
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parameter
        );

        void test();
        
        virtual void to_string();
        
        void setNurbsParameters(label rodNumber, label derivCoeffNumber, label dimension, scalar value);
        
        virtual vector rodDeriveParam(const LagrangianMarker* marker,const Parameter& par);
        virtual vector rodDeriveParam(label rodNumber,scalar rodParameter, scalar angle, scalar radiusFrac, const Parameter& par);

    protected:
        LineStructureStatus status;
        
        void initialize();
        
        virtual void check();
        
        virtual void createSpacingPoints();
        virtual void createSpacedPointsOnRod
        (
            label rodNumber,
            scalar spacing
        );

        void createMarkersFromSpacedPoints();
        virtual void createMarkersFromSpacedPointsOnRod
        (
            label rodNumber
        );

        void refineMarkers
        (
            std::pair<bool,scalar> forcedSpacing = {false,-1}
        );
        virtual void refineMarkersOnRod
        (
            label rodNumber,
            std::pair<bool,scalar> forcedSpacing = {false,-1}
        );

        void setMarkerVolume();
        virtual void setMarkerVolumeOnRod
        (
            label rodNumber
        );

        virtual void evaluateMarkerMeshRelation();
        void evaluateMarkerMeshRelation(std::list<LagrangianMarker>& markerList);

        virtual void reduceMarkers();

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
        );
        
        virtual void removeOverlapMarkers();

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
        
        virtual void collectMarkers();
        
        void collectHaloMarkers();
        
        void exchangeHaloMarkersData();
                
        void computeMarkerCellWeights();

        SystemSolve solutionStrategy = SystemSolve::Jacobi;
        using LinearSystem = std::tuple<CSR_Matrix_par,Vector_par>;
        std::unique_ptr<LinearSystem> computeMarkerEpsilonMatrix();
        void computeMarkerWeights();
        
        void exchangeHaloMarkersWeight();
        
                
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
        );

        template<typename T>
        static T integrateRodwise
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parA,
            scalar parB,
            std::function<T(scalar)> function
        );
        
        static scalar distance
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parA,
            scalar parB
        );
        
        static scalar distance(const LagrangianMarker& A, const LagrangianMarker& B);
        
        scalar bisectionBinary
        (
            scalar startValue,
            scalar endValue,
            std::function<bool(scalar)> criterion,
            scalar threshold
        );
        
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
    public:
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
                    std::tuple<vector,scalar,label,vector,DynamicList<Pair<label>>,DynamicList<vector>,DynamicList<scalar>,FixedList<scalar,10>,label,scalar,scalar,scalar> data
                );
                
                void insertMarkerWeight
                (
                    label haloCellInd,
                    label cellMarkerInd,
                    scalar weight
                );
                
                std::tuple<vector,scalar,label,vector,DynamicList<Pair<label>>,DynamicList<vector>,DynamicList<scalar>,FixedList<scalar,10>,label,scalar,scalar,scalar> getMarkerData
                (
                    label process,
                    label haloCellInd,
                    label cellMarkerInd
                ) const;
                
                scalar getMarkerWeight
                (
                    label process,
                    label haloCellInd,
                    label cellMarkerInd
                ) const;
                
                label size_processes() const;
                label size_haloCells(label process) const;
                label size_cellMarkers(label process, label haloCellInd) const;
                
                void communicate();
                
                void communicateWeight();
                
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
                List<List<DynamicList<FixedList<scalar,10>>>> globalHaloCellsMarkerb;
                // [nProc] -> [haloCells] -> [markers] -> weight
                List<List<DynamicList<scalar>>> globalHaloCellsMarkerWeight;
                // [nProc] -> [haloCells] -> [markers] -> nurbsInd
                List<List<DynamicList<label>>> globalHaloCellsMarkerNurbsInd;
                // [nProc] -> [haloCells] -> [markers] -> markerParameter
                List<List<DynamicList<scalar>>> globalHaloCellsMarkerParameter;
                // [nProc] -> [haloCells] -> [markers] -> markerAngle
                List<List<DynamicList<scalar>>> globalHaloCellsMarkerAngle;
                // [nProc] -> [haloCells] -> [markers] -> markerRadiusFrac
                List<List<DynamicList<scalar>>> globalHaloCellsMarkerRadiusFrac;
                
                void check();
                void check(label process, label haloCellInd, label cellMarkerInd) const;
                bool broadcasted;
                bool broadcastedWeights;
        };
    protected:
        // [haloCells] -> [marker of cell] -> ptr        
        std::vector<std::vector<std::pair<LagrangianMarker*,label>>> haloCellsRodMarkersList;
        //std::unique_ptr<List<List<DynamicList<scalar>>>> haloMarkerWeights;
        GlobalHaloMarkers globHaloMarkers;
    public:
        //inline std::unique_ptr<List<List<DynamicList<scalar>>>>& get_haloMarkerWeights(){return haloMarkerWeights;}
        inline const GlobalHaloMarkers& get_globalHaloMarkers(){return globHaloMarkers;}
        
        template<typename T>
        std::unique_ptr<List<List<DynamicList<T>>>> broadcastHaloMarkerField
        (
            const List<T>& collectedMarkerField
        );
        
        const std::vector<LagrangianMarker*>& getCollectedMarkers() const {return collectedMarkers;}
        
        //scalar initialSpacingFromMesh();
};
#include "FieldMarkerStructureInteraction.T.H"
#include "LineStructure.T.H"

class LineStructureParameters
{
    public:
        LineStructureParameters();
        virtual ~LineStructureParameters() = default;
        virtual void collectParameters(const LineStructure* structure);
        
    protected:
        std::vector<std::array<Parameter,3>> threeDimParameters;
};
}
#endif

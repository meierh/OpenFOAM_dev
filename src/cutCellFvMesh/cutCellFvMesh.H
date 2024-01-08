#ifndef CUTCELLFVMESH_H
#define CUTCELLFVMESH_H

//#include "fvCFD.H"
#include "dynamicRefineFvMesh.H"
#include <functional>
#include <memory>
#include "BsTree.H"
#include "KdTree.H"
#include "Nurbs1D.H"
#include "NurbsReader.H"
#include <unordered_set>
#include <unordered_map>
#include <list>
#include <algorithm>
#include <chrono>
#include <queue>
#include <cassert>
#include <bitset>
//#include <sys/types.h>
//#include <dirent.h>
//#include <string.h>
//#include <stdio.h>

#include "motionSolver.H"
#include "points0MotionSolver.H"
#include "displacementMotionSolver.H"
#include "valuePointPatchField.H"
#include "fixedValuePointPatchField.H"
#include "processorPointPatchField.H"
#include "motionDiffusivity.H"
#include "displacementLaplacianFvMotionSolver.H"
#include "patchWave.H"
#include "mapPolyMesh.H"
#include "polyTopoChange.H"

typedef double GRD_data_type;
typedef double MC33_real;

namespace Foam
{
//class motionSolver;
//class points0MotionSolver;
//class displacementMotionSolver;
    
class cutCellFvMesh   :   public dynamicRefineFvMesh //fvMesh
{
    //Public attributes
    public:
        enum cutStatus 
        {
            delNegMesh = 0,
            internalCut
        };
        
        struct nurbsReference
        {
            label nurbsInd=-1;
            scalar nurbsPara=-1;
        };
        
        
        const unsigned short int caseTable[2310] = 
        #include "MC33_LookUpTable.h"
        ;
        const std::array<std::array<unsigned short int,8>,24> permutationTable = {
        #include "CubePointPermutations.H"
        };
        const std::array<std::array<uint,2>,12> edgeToPnts = {{{0,1},{1,2},{2,3},{0,3},{4,5},{5,6},{6,7},{4,7},{0,4},{1,5},{2,6},{3,7}}};
        const std::unordered_map<uint,uint> pntsToEdge = {{1,0},{12,1},{23,2},{3,3},{45,4},{56,5},{67,6},{47,7},{4,8},{15,9},{26,10},{37,11}};
        const std::array<std::array<uint,4>,6> faceToPnts = {{{0,1,4,5},{1,2,5,6},{2,3,6,7},{0,3,4,7},{0,1,2,3},{4,5,6,7}}};
        const std::unordered_map<std::uint32_t,std::uint32_t> pntsToFace = {{145,0},{1256,1},{2367,2},{347,3},{123,4},{4567,5}};
        const std::unordered_map<std::uint8_t,std::unordered_set<std::uint8_t>> faceToEdges = {{0,{0,1,4,5}},{1,{1,5,9,10}},{2,{2,6,10,11}},{3,{3,7,8,11}},{4,{0,1,2,3}},{5,{4,5,6,7}}};
        
        void checkPermutation();

        class MC33{
            public:
                enum Case { c0=0,c1=1,c2=2,c31=3,c32=4,c411=5,c412=6,c5=7,c611=8,c612=9,c62=10,
                            c71=11,c72_720=12,c72_760=13,c72_800=14,c73_840=15,c73_912=16,c73_984=17,
                            c73=18,c741=19,c742=20,c8=21,c9=22,c1011_1190=23,c1011_1202=24,c1012_1214=25,
                            c1012_1238=26,c102_1262=27,c102_1286=28,c11=29,c1211_1358=30,c1211_1406=31,
                            c1212_1454=32,c1212_1550=33,c122_1646=34,c122_1742=35,c131=36,c132_MostNeg=37,
                            c132_MostPos=38,c133_MostPos=39,c133_MostNeg=40,c134=41,c1351=42,c1352_2206=43,
                            c1352_2246=45,c14=46
                            };
                enum permutationMethod{Point,Face,Both,FaceInversion};
                explicit MC33(cutCellFvMesh& mesh);
                
                typedef struct{
                    label cell=-1;
                    List<label> vertices = List<label>(8,-1);
                    List<edge> edges = List<edge>(12,edge(-1,-1));
                    List<label> edgeGlobalInd = labelList(12,-1);
                    List<label> cutEdgeVerticeIndex = labelList(12,-1);
                    label centerPointInd = -1;
                    std::unordered_map<label,label> cutVerticeIndexToEdge;
                    List<label> origFaces = List<label>(6);
                    std::vector<std::tuple<std::uint8_t,std::uint8_t,std::uint8_t>> cutTriangles;
                    unsigned int bitPattern;
                    std::array<std::int8_t,6> facePattern;
                    Case cubeCase = c0;
                    label permutationTableIndex;
                    bool redMarkIsPlusSide; //Efficient implementation of Marching Cubes cases with topological guarantees Figure 2
                    permutationMethod permMeth;
                    std::array<unsigned short int,8> pointPermutation;
                    std::array<std::uint8_t,12> edgePermutation;
                    std::array<std::uint8_t,6> facePermutation;
                } MC33Cube;
                
                MC33Cube computeCutCell(int cellInd);
                std::string getBitPattern(const unsigned short int x);
                
                
            private:
                cutCellFvMesh& mesh;
                const unsigned short int* table;
                const std::array<std::array<unsigned short int,8>,24>& permutations;
                std::array<std::array<unsigned short int,6>,24> facePermutations;
                const std::array<std::array<uint,2>,12> edgeToPnts;
                const std::unordered_map<uint,uint>& pntsToEdge;
                const std::array<std::array<uint,4>,6>& faceToPnts;
                const std::unordered_map<std::uint32_t,std::uint32_t>& pntsToFace;
                const std::unordered_map<std::uint8_t,std::unordered_set<std::uint8_t>>& faceToEdges;
                //const unsigned short int table[2310];
                
                std::array<std::int8_t,6> computeFacePattern();

                MC33_real v[8];
                //Procedures
                int face_tests(int *, int) const;
                int face_test1(int) const;
                int interior_test(int, int) const;

                unsigned int permuteBitPattern(const std::array<unsigned short int,8>& permutation,unsigned int bitPattern);
                void setBit(unsigned short int & bitField, std::uint8_t bitIndex);
                bool  getBit(const unsigned short int  bitField, std::uint8_t bitIndex);
                const unsigned short int* getTriangleCase(unsigned int verticePattern);
                std::vector<std::tuple<std::uint8_t,std::uint8_t,std::uint8_t>> collectTriangles
                (
                    const unsigned short int* triangleCase
                );
                MC33Cube generateMC33Cube(int cellInd);
                unsigned int computeSignBitPattern(const MC33Cube& cell);
                bool cellIsCut(int cellInd);
                bool detectInFaceTriangles(std::tuple<std::uint8_t,std::uint8_t,std::uint8_t>& triangle);
                void computeFacePermutations();
        };
        
        #include "ConvexCorrection.H"

    //Private Attributes
    private:
        // cellMinDimToStructureDim
        scalar cellDimToStructureDimLimit;
        
        // nurbs intersection radius
        scalar intersectionRadius = 0;
        
        //Motion Solver pointer
        autoPtr<motionSolver> motionPtr_;
        
        // Type of cut immersed boundary algorithm
        cutStatus ibAlgorithm;
        
        // Threshold of cell size difference
        scalar partialThreeshold = 0.1;
        
        // Level-Set Function for Boundary Computation
        std::function<scalar(const vector)> levelSet;
        
        // distance to Surface for all Points
        scalarList pointDist;
        
        // point vector to dist to reduce complexity
        std::unordered_map<label,label> newToOldPointIndMap;
        scalarList pointDistMap;
        
        
        //old Cell Volume 
        scalarList oldCellVolume;
        
        //new Cell Volume
        scalarList newCellVolume;
                
        //Points to Side (+1,0,-1)
        DynamicList<label> pointsToSide_;
        
        //Edge to Side (+1,0,-1)
        DynamicList<label> edgesToSide_;
        
        //face to Side (+1,0,-1)
        DynamicList<label> facesToSide_;
        
        //cell to Side (+1,0,-1)
        labelList cellsToSide_;
        
        //Set of all points on cut cell patch
        std::unordered_set<label> cutPatchPoint;
        
        //cutfaces to Side(+1,0,-1)
        DynamicList<label> cutFacesToSide_;
        
      //Points and data added by newMeshPoints function
      //Begin
            MC33 marchingCubesAlgorithm;
        
            //Cut data for every cell
            DynamicList<MC33::MC33Cube> mc33CutCellData;

            //Nbr of previous existing points
            label nbrOfPrevPoints;
            // new PointField; previous points + added Points
            DynamicList<point> newMeshPointsInFunc;
            pointField newMeshPoints_;
            List<bool> removedOldPoints;
        
            //Reference nurbs curves for zero points
            DynamicList<DynamicList<nurbsReference>> meshPointNurbsReference;
            DynamicList<point> meshPoint;
            
            //label of edge for every added Point and -1 for intersection point
            //labelList pointToEgde_;
            DynamicList<label> pointToEgde_;
            //label of added Point for every edge or -1 if edge has no added point
            //labelList edgeToPoint_;
            DynamicList<label> edgeToPoint_;
            
            //label of face for every added Point or zero side point and null for intersection non-cut point
            //labelListList pointToFaces_;
            //DynamicList<DynamicList<label>> pointToFaces_;
            //label of added points for every face or null if no point added
            //labelListList faceToPoints_;
            //DynamicList<DynamicList<label>> faceToPoints_;
            
            //label of cells for every added point and null for every intersection point
            //labelListList pointToCells_;
            //DynamicList<DynamicList<label>> pointToCells_;
            //label of added point for every cell or null if no added Points
            //labelListList cellToPoints_;
            //DynamicList<DynamicList<label>> cellToPoints_;
      //End
            
      //Edges and data added by addMeshEdges function
            //Begin
            //Nbr of previous existing edges
            label nbrOfPrevEdges;
            // new PointField; previous points + added Points
            DynamicList<edge> newMeshEdges_;
            List<bool> removedOldEdges;
        
            //label of face for every added edge or zero side edge and null every other edge
            DynamicList<DynamicList<label>> edgeToFaces_;
            //labelListList edgeToFaces_;
            //label of added points for every face or null if no point added
            DynamicList<DynamicList<label>> faceToEdges_;
            //labelListList faceToEdges_;
        
            //label of cells for every added point and null for every intersection point
            //DynamicList<DynamicList<label>> edgeToCells_;
            //labelListList edgeToCells_;
            //label of added edge or zero side edge for every cell or zero size
            //DynamicList<DynamicList<label>> cellToEdges_;
            //labelListList cellToEdges_;
            
            //true if face is problematic due to zero points -> 4,3 or 2 zero points
            //List<bool> problematicFace;
            //number of total zero points in problematic face
            //List<label> problematicFacePoints;
            //number of total new zero points in problematic face
            //List<label> problematicFaceNewPoints;
      //End
            
      //Faces and data added by addMeshFaces function
            //Begin
            //Nbr of previous existing edges
            label nbrOfPrevFaces;
            // new PointField; previous points + added Points
            DynamicList<face> newMeshFaces_;
            //faceList newMeshFaces_;
        
            //label of cells for every added face or zero side edge and null for every intersection point
            DynamicList<DynamicList<label>> faceToCells_;
            //labelListList faceToCells_;
            //label of added point for every cell or null if no added Points
            //DynamicList<DynamicList<label>> cellToFaces_;
            DynamicList<DynamicList<label>> cellToFaces_;
      //End
            
      //Split faces and data added by cutOldFaces
            //Begin
            // list of all cut faces replacing old faces
            DynamicList<face>    cutFaces_;
            // list of all faces containing a list of cut Faces or empty
            DynamicList<DynamicList<label>>   oldFacesToCutFaces_;
      //End
            
      //Data to construct new mesh
            //Begin
            // starting face for each patch
            labelList patchStarts;
            // nbr of faces for each patch
            labelList   patchSizes;
            
            //insert cell split faces
            // added cell split faces
            DynamicList<face>  addedCutFaces;
            // added split faces owner
            DynamicList<label> addedCutFacesOwner;
            // added split faces neighbor
            DynamicList<label> addedCutFacesNeighbor;
            // added split faces patch Index
            DynamicList<label> addedCutFacesPatchInd;
            
            //insert interior splitted faces
            // added interior split faces
            DynamicList<face>  splitAndUnsplitFacesInterior;
            // added interior split faces owner
            DynamicList<label> splitAndUnsplitFacesInteriorOwner;
            // added interior split faces neighbor
            DynamicList<label> splitAndUnsplitFacesInteriorNeighbor;
            // added interior split faces patch Index
            DynamicList<label> splitAndUnsplitFacesInteriorPatchInd;
            // added interior split faces old Index
            DynamicList<label> splitAndUnsplitFacesInteriorOldInd;
            DynamicList<label> splitAndUnsplitFacesInteriorOriginInd;
            
            //insert interior to boundary splitted faces
            // added interior to boundary split faces
            DynamicList<face>  splitAndUnsplitFacesInteriorToBoundary;
            // added interior to boundary split faces owner
            DynamicList<label> splitAndUnsplitFacesInteriorToBoundaryOwner;
            // added interior to boundary split faces neighbor
            DynamicList<label> splitAndUnsplitFacesInteriorToBoundaryNeighbor;
            // added interior to boundary split faces patch Index
            DynamicList<label> splitAndUnsplitFacesInteriorToBoundaryPatchInd;
            // added interior to boundary split faces old Index
            DynamicList<label> splitAndUnsplitFacesInteriorToBoundaryOldInd;
            DynamicList<label> splitAndUnsplitFacesInteriorToBoundaryOriginInd;
            
            //insert boundary splitted faces
            // added boundary split faces
            DynamicList<face>  splitAndUnsplitFacesBoundary;
            // added boundary split faces owner
            DynamicList<label> splitAndUnsplitFacesBoundaryOwner;
            // added boundary split faces neighbor
            DynamicList<label> splitAndUnsplitFacesBoundaryNeighbor;
            // added boundary split faces patch Index
            DynamicList<label> splitAndUnsplitFacesBoundaryPatchInd;
            // added boundary split faces old index
            DynamicList<label> splitAndUnsplitFacesBoundaryOldInd;
            DynamicList<label> splitAndUnsplitFacesBoundaryOriginInd;
            
            //insert interface splitted faces
            // added interface split faces
            DynamicList<face>  splitAndUnsplitFacesInterface;
            // added interface split faces owner
            DynamicList<label> splitAndUnsplitFacesInterfaceOwner;
            // added interface split faces neighbor
            DynamicList<label> splitAndUnsplitFacesInterfaceNeighbor;
            // added interface split faces patch Index
            DynamicList<label> splitAndUnsplitFacesInterfacePatchInd;
            // added interface split faces old Ind
            DynamicList<label>  splitAndUnsplitFacesInterfaceOldInd;
            DynamicList<label>  splitAndUnsplitFacesInterfaceOriginInd;
            
            //old face index to new face index
            std::unordered_map<label,DynamicList<label>> oldFaceToNewFaceInds;
            std::unordered_map<label,label> newFaceToOldFaceInds;
            
            //list of indexes of minus side cell of splitted cell or empty if cell is not split
            List<DynamicList<label>> oldSplittedCellToNewMinusCell;
            //list of indexes of plus side cell of splitted cell or empty if cell is not split
            List<DynamicList<label>> oldSplittedCellToNewPlusCell;
            //list of all cell, true for the ones that have been deleted
            List<bool> deletedCell;
            //list of labels to map old cell index to new plus side cut cell index; deleted Cells are set to -1
            DynamicList<DynamicList<label>> mapOldCellsToNewCells;
            //list of labels to map new cell index to old cells
            //labelList mapNewCellsToOldCells;
            
            //mapPolyMesh data
            label nOldPoints;
            labelList pointMap; 
            labelList purePointMap;
            labelList reversePointMap; 
            List<objectMap> pointsFromPoints; // not written
            pointField preMotionPoints;
            labelListList pointZoneMap; // not written
            
            label nOldFaces;
            labelList faceMap;
            labelList reverseFaceMap;
            List<objectMap> facesFromPoints; // not written
            List<objectMap> facesFromEdges; // not written
            List<objectMap> facesFromFaces; // not written
            labelListList faceZonePointMap; // not written
            labelListList faceZoneFaceMap; // not written
            labelHashSet flipFaceFlux; // not written
            
            label nOldCells;
            labelList cellMap;
            labelList reverseCellMap;
            List<objectMap> cellsFromPoints; // not written
            List<objectMap> cellsFromEdges; // not written
            List<objectMap> cellsFromFaces; // not written
            List<objectMap> cellsFromCells; // not written
            autoPtr<scalarField> oldCellVolumesPtr;
            labelListList cellZoneMap; // not written
            
            labelListList patchPointMap;
            labelList oldPatchStarts;
            labelList oldPatchNMeshPoints;
            
            std::unique_ptr<mapPolyMesh> meshCuttingMap;
        //End
            
        //Data for solidFraction Computation
            //index of faces for all cut cell faces
            labelListList cutCellsMinusAndPlus;
            // uncut cell to minus side cut cell face collection
            labelList oldCellToMinusCutCell;
            //uncut cell to plus side cut cell face collection
            labelList oldCellToPlusCutCell;
        //End
            
            labelList oldToNewPointInd;
            
            //New mesh data
            pointField new_points;
            faceList new_faces;
            labelList new_owner;
            labelList new_neighbour;
            
        //List<std::shared_ptr<Nurbs>> Curves;
        std::shared_ptr<std::vector<Nurbs1D>> Curves;
        //std::shared_ptr<std::vector<Nurbs1D>> Deformations;

        
        std::unique_ptr<KdTree> MainTree;
        List<std::unique_ptr<BsTree>> NurbsTrees;
        
        //std::unordered_map<label,scalar> pointToNurbsDist;
        
        
         
    //Public Memberfunctions    
    public:
        //Constructor
        explicit cutCellFvMesh
        (
            const IOobject& io,
            std::function<scalar(const vector)> levelSet,
            cutStatus state = delNegMesh
        );
        
        explicit cutCellFvMesh
        (
            const IOobject& io,
            std::shared_ptr<std::vector<Nurbs1D>> Curves,
            cutStatus state = delNegMesh
        );
        
        explicit cutCellFvMesh
        (
            const IOobject& io,
            std::shared_ptr<std::vector<Nurbs1D>> Curves,
            Time& runTime,
            std::unique_ptr<volScalarField>& solidFraction
        );
        
        explicit cutCellFvMesh
        (
            const IOobject& io,
            cutStatus state = delNegMesh,
            scalar cellDimToStructureDimLimit = 0.51
        );
        
        explicit cutCellFvMesh
        (
            const IOobject& io,
            scalar cellDimToStructureDimLimit = 0.51
        );

        explicit cutCellFvMesh
        (
            const IOobject& io,
            cutStatus state,
            NurbsReader& Reader,
            scalar cellDimToStructureDimLimit
        );
        
        void cutTheImmersedBoundary();
        void cutTheImmersedBoundary_MC33();
        void refineTheImmersedBoundary();
        void computeSolidFraction_MC33(std::unique_ptr<volScalarField>& solidFraction);

        
        // Print Functions for debugg purposes
        void printAddedPoints();
        void printAddedEdges();
        void printAddedFaces();
        void printCutFaces();
        void printNewMeshData();
        void printMesh();
        void selfTestMesh();
        
        //- Return the motionSolver
        const motionSolver& motion() const;

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
        
        // Set Defomration Curve
        void setInitialDeformationCurve
        (
            std::vector<scalarList>& nurbs_to_knots,
            std::vector<List<vector>>& nurbs_to_controlPoints,
            std::vector<scalarList>& nurbs_to_weights,
            std::vector<label>& nurbs_to_degree
        );
        
        //Compute immersed boundary mesh points movement vector
        void moveTheMesh();
        
        // move the Nurbs
        void moveNurbsCurves(std::vector<List<List<vector>>>& movedDeformationControlPoints);

        //- Write the underlying polyMesh and other data
        virtual bool writeObject
        (
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp,
            const bool write = true
        ) const;
        
        const DynamicList<DynamicList<nurbsReference>>& getMeshPointNurbsReference() const;
        
        const std::shared_ptr<std::vector<Nurbs1D>> getCurves() const;

    //Private Memberfunctions
    private:
        //Create pointsToSide Data
        void pointsToSide();
        void pointsToSide(const pointField& points);
        void findCutCellPnts();
        
        //Create edgesToSide Data
        void edgesToSide();
        void edgesToSide(const edgeList& edges);
        
        //Create facesToSide Data
        void facesToSide();
        void facesToSide(const faceList& faces);
        
        //Create cellsToSide Data
        void cellsToSide();
        void cellsToSide(const cellList& cells);
        
        //Create Marching Cube data
        void executeMarchingCubes();
        
        // Add points for cutting the cells and create dependency data for creating cut cells
        // Adds to points_.
        void newMeshPoints();
        void newMeshPoints_MC33();
        
        scalar distToNurbs(point pnt,bool& foundFlag,DynamicList<nurbsReference>& reference);
        scalar distToNurbs(point pnt,bool& foundFlag);

        label sideToNurbs(point pnt,bool& foundFlag);
        
        //std::unique_ptr<mapPolyMesh> CutCellPolyMeshMap;
        
        // Add edges for cutting the cells and create dependency data for creating cut cells
        // Adds to edges_;
        void newMeshEdges_MC33();
        
        
        bool equalEdges
        (
            std::unordered_set<label> setA,
            std::unordered_set<label> setB
        );
        /*
        void findCycles
        (
            int i,
            label startPoint,
            label nextPoint,
            label prevPoint,
            DynamicList<label> cyclePath,
            DynamicList<label> cycleEdgePath,
            std::unordered_set<label> coveredPoints,
            std::unordered_set<label> usedEdges,
            std::unordered_map<label,std::unordered_set<label>>& pointGraphData,
            DynamicList<std::pair<label,label>>& pointEdgeComb,
            DynamicList<DynamicList<label>>& closedCyclePoints,
            DynamicList<std::unordered_set<label>>& closedCycleEdges,
            DynamicList<DynamicList<label>>& closedCycleEdgesList
        );
        */
        vector vectorToNurbs
        (
            point pnt,
            bool& foundFlag
        );
        List<vector> vectorsToNurbsOfEdge
        (
            point startPoint,
            point endPoint,
            label nbrOfVectors = 10
        );
        
        List<scalar> distToNursOfEdge
        (
            point startPoint,
            point endPoint,
            label nbrOfPoints = 10
        );

        scalar nearestNurbsIndexPara
        (
            point pnt,
            bool& foundFlag,
            label& nurbsInd,
            scalar& nurbsPara
        );
        
        /*
        void computeClosedFaceFront
        (
            label centerPointInd,
            label faceInd,
            labelList& problematicFacePoints,
            const labelList& cellsInd,
            List<DynamicList<face>>& facesInCellsIn,
            List<DynamicList<std::unordered_set<label>>>& facesMapInCellsIn,
            DynamicList<bool>& closedFaceFrontOutFaceInFace,
            DynamicList<DynamicList<face>>& closedFaceFrontOut,
            DynamicList<DynamicList<std::unordered_set<label>>>& closedFaceFrontMapOut
        );
        */
        
        /*Greg Lukianoff
        void groupEdgesToFaces
        (
            List<DynamicList<DynamicList<label>>>& cellNonconnectedEdges,
            List<DynamicList<DynamicList<DynamicList<label>>>>& cellNonConnectedMultiPoints,
            List<DynamicList<DynamicList<DynamicList<label>>>>& cellNonConnectedMultiEdges,
            List<DynamicList<DynamicList<face>>>& cellNonConnectedMultiFaces,
            List<DynamicList<DynamicList<std::unordered_set<label>>>>& cellNonConnectedMultiPointMap,
            List<DynamicList<DynamicList<std::unordered_set<label>>>>& cellNonConnectedMultiEdgeMap
        );
        */
        
        // Add faces for cutting the cells and create dependency data for creating cut cells
        // Adds to faces_;
        void newMeshFaces_MC33();
        
        // Splits all old faces being cut into too
        void cutOldFaces_MC33();

        //rebuild mesh deleting negative side
        void createNewMeshData_MC33();
        
        // project distance from Nurbs Surface on Points
        void projectNurbsSurface(bool reset=false);
        
        // project level set function on Points
        void projectLevelSet();
        
        // correct the point order in faces to correct the normal dir 
        void correctFaceNormalDir(const pointField& points,faceList& faces,const labelList& owner,const labelList& neighbour);
        List<bool> correctListFaceNormalDir(const pointField& points,const faceList& faces,const labelList& owner,const labelList& neighbour);

        //cell agglomeration
        /*
        void agglomerateSmallCells_cutNeg
        (
            scalarList& newCellVolume,
            scalarList& oldCellVolume,
            scalar partialThreeshold = 0.25
        );
        */
        
        typedef struct{
            labelList originalFaceInds;
            labelList splittedFaceInds;
            labelList replacingFaceInds;
            labelList replacingResidualFaceInds;
            labelList addedFaceInds;
        } CellFaces;
        
        typedef struct{
            label originalCell = -1;
            DynamicList<std::pair<vector,edge>> addedPoints;
            label addedPointsLimit=-1;
            // face with possibly added Points as index+addedPointsLimit
            DynamicList<face> addedFace;
            // face with reference index for original face
            DynamicList<std::pair<face,label>> splittedFaces;
            // face with reference index for replacing face and original(replaced) face pair(orig or splitted = true / index)
            DynamicList<std::tuple<face,label,std::pair<bool,label>,bool>> replacingFaces;
            // face with reference index for replacing face
            DynamicList<std::pair<face,label>> replacingResidualFaces;
            // new cells with tuple of (original face indexes, splitted faces indexex, added faces indexes)
            DynamicList<CellFaces> cells;
            
            DynamicList<label> nonCellFaceChanges;
        } CellSplitData;
        
        label memSecMC33VerticePntInd
        (
            label mc33LocalPntInd_I,
            const MC33::MC33Cube& mc33cube,
            const std::unordered_set<label>& verticeOfCell
        );
        
        label memSecMC33CutEdgePntInd
        (
            label mc33LocalEdgeInd_I,
            const MC33::MC33Cube& mc33cube,
            const std::unordered_set<label>& verticeOfCell
        );
        
        label getFaceIndFromPntList
        (
            const cell& thisCell,
            std::vector<label> pointList
        );
        
        void triangulateCell
        (
            const cell& thisCell
        );
        
        std::unique_ptr<List<scalar>> faceToCellCenterRelation
        (
            const cell& thisCell,
            const label cellInd
        );
        
        std::unique_ptr<List<DynamicList<std::pair<label,std::pair<label,label>>>>> compCellFaceEdgeGraph
        (
            const cell& oneCell                                                                             
        );
        std::unique_ptr<List<List<DynamicList<label>>>> compCellFacePointGraph
        (
            const cell& oneCell                                                                                 
        );
        
        bool testCorrectFaceSizes
        (
            const std::array<label,10> faceSizeNbrReq,
            const cell& thisCell,
            bool print=false
        );
        
        std::unique_ptr<cellList> createCells
        (
            const faceList& faces,
            const labelList& owner,
            const labelList& neighbour 
        );
        
        void testCellClosure
        (
            const faceList& faces,
            const cellList& cells
        );
        
        void correctFaceNormal
        (
            faceList& faces,
            const labelList& owner,
            const labelList& neighbour,
            const cellList& cells
        );
        
        void agglomerateSmallCells_MC33
        (
            scalar partialThreeshold = 0.01
        );
        
        bool verticesAreIdentical
        (
            const corrFaceVertice& vert1,
            const corrFaceVertice& vert2,
            const MC33::MC33Cube& mc33cube,
            corrFaceVertice& vertResult
        );
        
        void getGlobalPoint
        (
            const corrFaceVertice& vert,
            const MC33::MC33Cube& mc33cube,
            const pointField& points,
            label& verticeInd,
            vector& verticePoint
        );

        label getGlobalPoint
        (
            const corrFaceVertice& vert,
            const MC33::MC33Cube& mc33cube,
            const std::unordered_map<std::uint8_t,label>& edgeToAddedPntInd,
            const pointField& points
        );

        face abstrFaceToGlobalFaceTransfer
        (
            const corrFace& abstrFace,
            const MC33::MC33Cube& mc33cube,
            const std::unordered_map<std::uint8_t,label>& edgeToAddedPntInd,
            const faceList& faces,
            const pointField& points
        );
        
        enum FaceRelation{Ident=0,oneInTwo=1,twoInOne=2,mixedOverlay=3,none=4};
        FaceRelation facesOverlay
        (
            const face& face1,
            const face& face2
        );
        
        std::unique_ptr<CellSplitData> generateCellSplitData
        (
            const cell& thisCell,
            const label cellInd,
            corrConcaveCellData& corrData,
            std::function<bool(const corrConcaveCellData&, const List<face>&, label)> isCellValid =
            [] (const corrConcaveCellData& corrData, const List<face>& nonConvexCellFaces, label cellInd)
            {
                const List<label>& oneCell = corrData.cells[cellInd];
                List<corrFace> oneCellAbstrFaces(oneCell.size());
                std::transform(oneCell.begin(),oneCell.end(),oneCellAbstrFaces.begin(),[&](label faceInd){return corrData.faces[faceInd];});
                List<face> oneCellGlobFaces(oneCell.size());
                std::transform(oneCell.begin(),oneCell.end(),oneCellGlobFaces.begin(),[&](label faceInd){return nonConvexCellFaces[faceInd];});
                
                List<DynamicList<label>> faceInFace(oneCellGlobFaces.size());
                List<DynamicList<label>> faceIdent(oneCellGlobFaces.size());
                for(uint baseInd=0; baseInd<oneCellGlobFaces.size(); baseInd++)
                {
                    face baseFace = oneCellGlobFaces[baseInd];
                    std::unordered_set<label> baseFaceSet(baseFace.begin(),baseFace.end());
                    for(uint compareInd=0; compareInd<oneCellGlobFaces.size(); compareInd++)
                    {
                        face compFace = oneCellGlobFaces[compareInd];
                        std::unordered_set<label> compFaceSet(compFace.begin(),compFace.end());
                        if(baseInd==compareInd)
                            continue;
                        bool allIn=true;
                        for(label compVert : compFace)
                        {
                            if(baseFaceSet.find(compVert)==baseFaceSet.end())
                                allIn=false;
                        }
                        if(allIn)
                        {
                            if(baseFace.size()==compFace.size())
                                faceIdent[baseInd].append(compareInd);
                            else if(baseFace.size()==compFace.size())
                                faceInFace[baseInd].append(compareInd);
                        }
                    }
                }
                std::unordered_set<label> facesInFaceSet;
                for(DynamicList<label> inFaceList : faceInFace)
                {
                    for(label inFace : inFaceList)
                    {
                        if(facesInFaceSet.find(inFace)!=facesInFaceSet.end())
                            FatalErrorInFunction<<"Face can not be in two base faces!"<< exit(FatalError);
                        facesInFaceSet.insert(inFace);
                    }
                }
                
                DynamicList<face> reducedOneCell;
                for(label ind=0; ind<oneCellGlobFaces.size(); ind++)
                {
                    face oneFace = oneCellGlobFaces[ind];
                    if(facesInFaceSet.find(ind)==facesInFaceSet.end())
                    {
                        DynamicList<label>& ident = faceIdent[ind];
                        if(ident.size()>0)
                        {
                            std::sort(ident.begin(),ident.end());
                            if(ind<ident[0])
                                reducedOneCell.append(oneFace);
                        }
                        else
                            reducedOneCell.append(oneFace);
                    }
                }
                                    
                std::unordered_set<label> uniquePoints;
                for(face& oneGlobFace : reducedOneCell)
                    for(label globVert : oneGlobFace)
                        uniquePoints.insert(globVert);
                
                if(reducedOneCell.size()<4)
                    return false;
                List<label> reducedOneCellFaceSize(reducedOneCell.size());
                std::transform(reducedOneCell.begin(),reducedOneCell.end(),reducedOneCellFaceSize.begin(),
                               [&](const face& oneFace){return oneFace.size();});
                std::sort(reducedOneCellFaceSize.begin(),reducedOneCellFaceSize.end(),std::greater<label>());
                
                if(reducedOneCell.size() < reducedOneCellFaceSize[0]+1)
                    return false;
                return true;
            }
        );
        
        void testCellSplitData
        (
            const cell& thisCell,
            const faceList& faces,
            const pointField& points,
            const CellSplitData& data
        );
        
        bool redMarkSide(const MC33::MC33Cube& mc33cube);
        
        std::unique_ptr<CellSplitData> nonConvexCellSplitC2
        (
            const cell& thisCell,
            const label cellInd,
            std::unique_ptr<List<scalar>> faceToCellCenterRelation
        );
        
        std::unique_ptr<CellSplitData> nonConvexCellSplitC32
        (
            const cell& thisCell,
            const label cellInd,
            std::unique_ptr<List<scalar>> faceToCellCenterRelation
        );
        
        std::unique_ptr<CellSplitData> nonConvexCellSplitC412
        (
            const cell& thisCell,
            const label cellInd,
            std::unique_ptr<List<scalar>> faceToCellCenterRelation
        );
        
        std::unique_ptr<CellSplitData> nonConvexCellSplitC5
        (
            const cell& thisCell,
            const label cellInd,
            std::unique_ptr<List<scalar>> faceToCellCenterRelation
        );
        
        std::unique_ptr<CellSplitData> nonConvexCellSplitC612
        (
            const cell& thisCell,
            const label cellInd,
            std::unique_ptr<List<scalar>> faceToCellCenterRelation
        );
        
        std::unique_ptr<CellSplitData> nonConvexCellSplitC62
        (
            const cell& thisCell,
            const label cellInd,
            std::unique_ptr<List<scalar>> faceToCellCenterRelation
        );
        
        std::unique_ptr<CellSplitData> nonConvexCellSplitC72
        (
            const cell& thisCell,
            const label cellInd,
            std::unique_ptr<List<scalar>> faceToCellCenterRelation
        );
        
        std::unique_ptr<CellSplitData> nonConvexCellSplitC73
        (
            const cell& thisCell,
            const label cellInd,
            std::unique_ptr<List<scalar>> faceToCellCenterRelation
        );
        
        std::unique_ptr<CellSplitData> nonConvexCellSplitC8
        (
            const cell& thisCell,
            const label cellInd,
            std::unique_ptr<List<scalar>> faceToCellCenterRelation
        );
        
        std::unique_ptr<CellSplitData> nonConvexCellSplitC1011
        (
            const cell& thisCell,
            const label cellInd,
            std::unique_ptr<List<scalar>> faceToCellCenterRelation
        );
        
        std::unique_ptr<CellSplitData> nonConvexCellSplitC1211
        (
            const cell& thisCell,
            const label cellInd,
            std::unique_ptr<List<scalar>> faceToCellCenterRelation
        );
        
        std::unique_ptr<CellSplitData> nonConvexCellSplitC122
        (
            const cell& thisCell,
            const label cellInd,
            std::unique_ptr<List<scalar>> faceToCellCenterRelation
        );
        
        std::unique_ptr<CellSplitData> nonConvexCellSplitC132
        (
            const cell& thisCell,
            const label cellInd,
            std::unique_ptr<List<scalar>> faceToCellCenterRelation
        );
        
        std::unique_ptr<CellSplitData> nonConvexCellSplitC133
        (
            const cell& thisCell,
            const label cellInd,
            std::unique_ptr<List<scalar>> faceToCellCenterRelation
        );
        
        std::unique_ptr<CellSplitData> nonConvexCellSplitC134
        (
            const cell& thisCell,
            const label cellInd,
            std::unique_ptr<List<scalar>> faceToCellCenterRelation
        );
        
       std::unique_ptr<CellSplitData> nonConvexCellSplitC14
        (
            const cell& thisCell,
            const label cellInd,
            std::unique_ptr<List<scalar>> faceToCellCenterRelation
        );
        
        void splitFaceByEdge
        (
            const face& oneFace,
            const edge& oneEdge,
            FixedList<face,2>& splitFaces
        );
        
        bool edgeInFace
        (
            const face& oneFace,
            const edge& oneEdge
        );
        
        void facesEdgesIntersection
        (
            const face& totalFace,
            const face& ownFace,
            const face& neiFace,
            List<List<bool>>& faceEdgeIntersections
        );
        
        bool faceEdgesIntersect
        (
            const List<List<bool>>& faceEdgeIntersections
        );
        
        vector computeIntersectionPnt
        (
            const edge edgeA,
            const edge edgeB
        );
        
        struct ArrayHasher {
            std::size_t operator()(const FixedList<label,4>& a) const {
                std::size_t h = 0;
                for (label e : a) {
                    h ^= std::hash<label>{}(e)  + 0x9e3779b9 + (h << 6) + (h >> 2); 
                }
                return h;
            }
        };
        
        // (pnt0,pnt1),(pnt2,pnt3) as key with (pnt0<pnt1) and (pnt2<pnt3) and pnt0<pnt2
        using edgToIntPntIndMap = std::unordered_map<FixedList<label,4>,label,ArrayHasher>;
        List<label> faceIntersection
        (
            const face& totalFace,
            const face& ownFace,
            const face& neiFace,
            edgToIntPntIndMap& edgesToAddPntInd,
            DynamicList<vector>& addPnt
        );
        
        edge getEdge
        (
            const face& theFace,
            label edgeInd
        );
        
        std::unique_ptr<std::vector<std::pair<label,label>>> matchingPoints
        (
            const face& ownFace,
            const face& neiFace
        );

        std::array<std::pair<label,label>,2> matchingEdge
        (
            const face& ownFace,
            const face& neiFace
        );
        
        scalar computeCellSizeAndSizeLimit
        (
            const List<cell>& cells,
            List<bool> tooSmallCells,
            List<scalar> cellSizes,
            scalar limitPercentage = 0.1,
            scalar tooSmallFraction = 0.1
        );
        
        void findSmallCellMergeCandidate
        (
            const List<cell>& cells,
            const List<bool> smallCell,
            List<DynamicList<std::pair<label,label>>>& smallCellMergeCand
        );

        void testCellConvexity
        (
            const List<cell>& cells,
            List<List<DynamicList<std::tuple<label,std::pair<label,label>,bool,scalar>>>>& cellFaceEdgeGraph,
            List<bool>& nonConvexCell
        );
        
        /*
        void agglomerateSmallCells_cutNeg_plus
        (
            scalarList& newCellVolume,
            scalarList& oldCellVolume,
            scalar partialThreeshold = 0.25
        );
        */
        
        /*
        labelList searchDown
        (
            scalarListList& possibleMergeFaceArea,
            labelListList& possibleMergeFaces,
            labelListList& possibleMergeCells,
            boolList& oneMergeFaceSufficient,
            boolList& mergeNecessary,
            label count,
            std::unordered_set<label> cellReserved
        );
        */
        
        /*
        labelList searchDown_rec
        (
            DynamicList<DynamicList<scalar>>& possibleMergeFaceArea,
            DynamicList<DynamicList<label>>& possibleMergeFaces,
            DynamicList<DynamicList<label>>& possibleMergeCells,
            DynamicList<bool>& oneMergeFaceSufficient,
            DynamicList<bool>& mergeNecessary,
            label count,
            DynamicList<DynamicList<label>>& blockedCells,
            std::unordered_set<label>& cellReserved
        );
        */
        
        /*
        List<DynamicList<label>> searchDown_iter
        (
            DynamicList<DynamicList<scalar>>& possibleMergeFaceArea,
            DynamicList<DynamicList<DynamicList<label>>>& possibleMergeFaces,
            DynamicList<DynamicList<DynamicList<label>>>& possibleMergeCells,
            DynamicList<bool>& oneMergeFaceSufficient,
            DynamicList<bool>& mergeNecessary
        );
        */
        
        /*
        List<DynamicList<label>> searchDown_iter_preBlock
        (
            const labelList& owner,
            const labelList& neighbour,
            DynamicList<DynamicList<scalar>>& possibleMergeFaceArea,
            DynamicList<DynamicList<DynamicList<label>>>& possibleMergeFaces,
            DynamicList<DynamicList<DynamicList<label>>>& possibleMergeCells,
            DynamicList<bool>& oneMergeFaceSufficient,
            DynamicList<bool>& mergeNecessary
        );
        */
        
        /*
        List<DynamicList<label>> assignMergeFaces
        (
            const labelList& owner,
            const labelList& neighbour,
            DynamicList<DynamicList<scalar>>& possibleMergeFaceArea,
            DynamicList<DynamicList<DynamicList<label>>>& possibleMergeFaces,
            DynamicList<DynamicList<DynamicList<label>>>& possibleMergeCells,
            DynamicList<bool>& oneMergeFaceSufficient,
            DynamicList<bool>& mergeNecessary,
            DynamicList<DynamicList<scalar>>& possibleMergeCellsPartialSize,
            scalar partialThreeshold
        );
        */
        
        //test data for new Mesh creation
        void testNewMeshData
        (
            const faceList& newFaces,
            const labelList& newFaceOwner,
            const labelList& newFaceNeighbor,
            const labelList& patchStarts,
            const labelList& patchSizes
        );
        
        /*
        void testForCellSize
        (
            DynamicList<DynamicList<scalar>>& possibleMergeFaceArea,
            DynamicList<DynamicList<DynamicList<label>>>& possibleMergeFaces,
            DynamicList<DynamicList<DynamicList<label>>>& possibleMergeCells,
            DynamicList<bool>& oneMergeFaceSufficient,
            DynamicList<bool>& mergeNecessary,
            labelList& mergeFaceOfCell,
            scalar partialThreeshold
        );
        */
        
        void checkForHexCellsInCutArea();
                
        scalar minNurbsRadius();
        
        void Barrier(bool stop=false);
};
void UnitTest_cutCellFvMesh(int argc, char *argv[], Time& runtime);
};

#endif
